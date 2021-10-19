# ======================================================================
#         Imports
# ======================================================================
from collections import OrderedDict
import numpy
from mpi4py import MPI
from baseclasses.utils import Error
from pysurf import intersectionAPI, curveSearchAPI, utilitiesAPI, adtAPI, tsurf_tools, tecplot_interface
from . import DVGeometry


class DVGeometryMulti:
    """
    A class for manipulating multiple components using multiple FFDs

    """

    def __init__(self, comm=MPI.COMM_WORLD, dh=1e-6):

        self.compNames = []
        self.comps = OrderedDict()
        self.DVGeoDict = OrderedDict()
        self.points = OrderedDict()
        self.comm = comm
        self.updated = {}
        self.dh = dh
        self.intersectComps = []

        # flag to keep track of IC jacobians
        self.ICJupdated = False

    def addComponent(self, comp, ffdFile, triMesh=None, scale=1.0, bbox={}):
        """
        Method to add components to the DVGeometryMulti object.
        Returns the DVGeo object for this component
        """

        # we need to create a new DVGeo object for this component
        DVGeo = DVGeometry(ffdFile)

        if triMesh is not None:
            # We also need to read the triMesh and save the points
            nodes, triConn, barsConn = self._readCGNSFile(triMesh)

            # scale the nodes
            nodes *= scale

            # add these points to the corresponding dvgeo
            DVGeo.addPointSet(nodes, "triMesh")
        else:
            # the user has not provided a triangulated surface mesh for this file
            nodes = None
            triConn = None
            barsConn = None

        # we will need the bounding box information later on, so save this here
        xMin, xMax = DVGeo.FFD.getBounds()

        # also we might want to modify the bounding box if the user specified any coordinates
        if "xmin" in bbox:
            xMin[0] = bbox["xmin"]
        if "ymin" in bbox:
            xMin[1] = bbox["ymin"]
        if "zmin" in bbox:
            xMin[2] = bbox["zmin"]
        if "xmax" in bbox:
            xMax[0] = bbox["xmax"]
        if "ymax" in bbox:
            xMax[1] = bbox["ymax"]
        if "zmax" in bbox:
            xMax[2] = bbox["zmax"]

        # initialize the component object
        self.comps[comp] = component(comp, DVGeo, nodes, triConn, barsConn, xMin, xMax)

        # add the name to the list
        self.compNames.append(comp)

        # also save the DVGeometry pointer in the dictionary we pass back
        self.DVGeoDict[comp] = DVGeo

        return DVGeo

    def addIntersection(
        self,
        compA,
        compB,
        dStarA=0.2,
        dStarB=0.2,
        featureCurves=[],
        distTol=1e-14,
        project=False,
        marchDir=1,
        includeCurves=False,
        intDir=None,
        curveEpsDict={},
        remeshBwd=True,
    ):
        """
        Method that defines intersections between components
        """

        # just initialize the intersection object
        self.intersectComps.append(
            CompIntersection(
                compA,
                compB,
                dStarA,
                dStarB,
                featureCurves,
                distTol,
                self,
                project,
                marchDir,
                includeCurves,
                intDir,
                curveEpsDict,
                remeshBwd,
            )
        )

    def getDVGeoDict(self):
        # return DVGeo objects so that users can add design variables
        return self.DVGeoDict

    def finalizeDVs(self):
        """
        This function should be called after adding all DVGeoDVs
        """

        self.DV_listGlobal = OrderedDict()  # Global Design Variable List
        self.DV_listLocal = OrderedDict()  # Local Design Variable List
        self.DV_listSectionLocal = OrderedDict()  # Local Normal Design Variable List

        # we loop over all components and add the dv objects
        for comp in self.compNames:

            # get this DVGeo
            DVGeoComp = self.comps[comp].DVGeo

            # loop over the DVGeo's DV lists

            for k, v in DVGeoComp.DV_listGlobal.items():
                # change the key and add it to our dictionary...
                knew = comp + ":" + k
                self.DV_listGlobal[knew] = v

            for k, v in DVGeoComp.DV_listLocal.items():
                # change the key and add it to our dictionary...
                knew = comp + ":" + k
                self.DV_listLocal[knew] = v

            for k, v in DVGeoComp.DV_listSectionLocal.items():
                # change the key and add it to our dictionary...
                knew = comp + ":" + k
                self.DV_listSectionLocal[knew] = v

    def addPointSet(self, points, ptName, compNames=None, comm=None, applyIC=False, **kwargs):

        # if the user passes a list of compNames, we only use these comps.
        # we will still use all points to add the pointset, but by default,
        # the components not in this list will get 0 npoints. This is for
        # consistency, each dvgeo will have all pointsets and book keeping
        # becomes easier this way.

        # if compList is not provided, we use all components
        if compNames is None:
            compNames = self.compNames

        # before we do anything, we need to create surface ADTs
        # for which the user provided triangulated meshes
        # TODO Time these, we can do them once and keep the ADTs
        for comp in compNames:

            # check if we have a trimesh for this component
            if self.comps[comp].triMesh:

                # now we build the ADT
                # from ney's code:
                # Set bounding box for new tree
                BBox = numpy.zeros((2, 3))
                useBBox = False

                # dummy connectivity data for quad elements since we have all tris
                quadConn = numpy.zeros((0, 4))

                # Compute set of nodal normals by taking the average normal of all
                # elements surrounding the node. This allows the meshing algorithms,
                # for instance, to march in an average direction near kinks.
                nodal_normals = adtAPI.adtapi.adtcomputenodalnormals(
                    self.comps[comp].nodes.T, self.comps[comp].triConn.T, quadConn.T
                )
                self.comps[comp].nodal_normals = nodal_normals.T

                # Create new tree (the tree itself is stored in Fortran level)
                adtAPI.adtapi.adtbuildsurfaceadt(
                    self.comps[comp].nodes.T,
                    self.comps[comp].triConn.T,
                    quadConn.T,
                    BBox.T,
                    useBBox,
                    MPI.COMM_SELF.py2f(),
                    comp,
                )

        # create the pointset class
        self.points[ptName] = PointSet(points, comm=comm)

        for comp in self.compNames:
            # initialize the list for this component
            self.points[ptName].compMap[comp] = []
            self.points[ptName].compMapFlat[comp] = []

        # we now need to create the component mapping information
        for i in range(self.points[ptName].nPts):

            # initial flags
            inFFD = False
            proj = False
            projList = []

            # loop over components and check if this point is in a single BBox
            for comp in compNames:

                # apply a small tolerance for the bounding box in case points are coincident with the FFD
                boundTol = 1e-16
                xMin = self.comps[comp].xMin
                xMax = self.comps[comp].xMax
                xMin -= numpy.abs(xMin * boundTol) + boundTol
                xMax += numpy.abs(xMax * boundTol) + boundTol

                # check if inside
                if (
                    xMin[0] < points[i, 0] < xMax[0]
                    and xMin[1] < points[i, 1] < xMax[1]
                    and xMin[2] < points[i, 2] < xMax[2]
                ):

                    # print('point',i,'is in comp',comp)
                    # add this component to the projection list
                    projList.append(comp)

                    # this point was not inside any other FFD before
                    if not inFFD:
                        inFFD = True
                        inComp = comp
                    # this point was inside another FFD, so we need to project it...
                    else:
                        # set the projection flag
                        proj = True

            # project this point to components, we need to set inComp string
            if proj:

                # set a high initial distance
                dMin2 = 1e10

                # loop over the components
                for comp in compNames:
                    # check if this component is in the projList
                    if comp in projList:

                        # check if we have an ADT:
                        if self.comps[comp].triMesh:
                            # Initialize reference values (see explanation above)
                            numPts = 1
                            dist2 = numpy.ones(numPts) * 1e10
                            xyzProj = numpy.zeros((numPts, 3))
                            normProjNotNorm = numpy.zeros((numPts, 3))

                            # Call projection function
                            _, _, _, _ = adtAPI.adtapi.adtmindistancesearch(
                                points[i].T, comp, dist2, xyzProj.T, self.comps[comp].nodal_normals.T, normProjNotNorm.T
                            )
                            # print('Distance of point',points[i],'to comp',comp,'is',numpy.sqrt(dist2))
                            # if this is closer than the previous min, take this comp
                            if dist2 < dMin2:
                                dMin2 = dist2[0]
                                inComp = comp

                        else:
                            raise Error(
                                "The point at \n(x, y, z) = (%.3f, %.3f, %.3f) \nin pointset %s is inside multiple FFDs but a triangulated mesh for component %s is not provided to determine which component owns this point."
                                % (points[i, 0], points[i, 1], points[i, 1], ptName, comp)
                            )

            # this point was inside at least one FFD. If it was inside multiple,
            # we projected it before to figure out which component it should belong to
            if inFFD:
                # we can add the point index to the list of points inComp owns
                self.points[ptName].compMap[inComp].append(i)

                # also create a flattened version of the compMap
                for j in range(3):
                    self.points[ptName].compMapFlat[inComp].append(3 * i + j)

            # this point is outside any FFD...
            else:
                raise Error(
                    "The point at (x, y, z) = (%.3f, %.3f, %.3f) in pointset %s is not inside any FFDs"
                    % (points[i, 0], points[i, 1], points[i, 1], ptName)
                )

        # using the mapping array, add the pointsets to respective DVGeo objects
        for comp in self.compNames:
            compMap = self.points[ptName].compMap[comp]
            # print(comp,compMap)
            self.comps[comp].DVGeo.addPointSet(points[compMap], ptName)

        # check if this pointset will get the IC treatment
        if applyIC:
            # loop over the intersections and add pointsets
            for IC in self.intersectComps:
                IC.addPointSet(points, ptName, self.points[ptName].compMap, comm)

        # finally, we can deallocate the ADTs
        for comp in compNames:
            if self.comps[comp].triMesh:
                adtAPI.adtapi.adtdeallocateadts(comp)
                # print('Deallocated ADT for component',comp)

        # mark this pointset as up to date
        self.updated[ptName] = False

    def setDesignVars(self, dvDict):
        """
        Standard routine for setting design variables from a design
        variable dictionary.

        Parameters
        ----------
        dvDict : dict
            Dictionary of design variables. The keys of the dictionary
            must correspond to the design variable names. Any
            additional keys in the dfvdictionary are simply ignored.
        """

        # first get the list of DVs from each comp so we can ignore extra entries
        for comp in self.compNames:
            self.comps[comp].dvDict = self.comps[comp].DVGeo.getValues()

        # loop over all dvs we get as the input
        for k, v in dvDict.items():
            # we only set dvgeomulti DVs. Then k should always have a : in it
            if ":" in k:
                # get the component name
                comp, dvName = k.split(":", 1)

                # now check if this comp has this dv
                if dvName in self.comps[comp].dvDict:
                    # set the value
                    self.comps[comp].dvDict[dvName] = v

        # loop over the components and set the values
        for comp in self.compNames:
            self.comps[comp].DVGeo.setDesignVars(self.comps[comp].dvDict)

        # We need to give the updated coordinates to each of the
        # intersectComps (if we have any) so they can update the new
        # intersection curve
        for IC in self.intersectComps:
            IC.setSurface(self.comm)

        # Flag all the pointSets as not being up to date:
        for pointSet in self.updated:
            self.updated[pointSet] = False

        # also set IC Jacobians as out of date
        self.ICJupdated = False

    def getValues(self):
        """
        Generic routine to return the current set of design
        variables. Values are returned in a dictionary format
        that would be suitable for a subsequent call to setValues()

        Returns
        -------
        dvDict : dict
            Dictionary of design variables
        """

        dvDict = {}
        # we need to loop over each DVGeo object and get the DVs
        for comp in self.compNames:
            dvDictComp = self.comps[comp].DVGeo.getValues()
            # we need to loop over these DVs.
            for k, v in dvDictComp.items():
                # We will add the name of the comp and a : to the full DV name
                dvName = "%s:%s" % (comp, k)
                dvDict[dvName] = v

        return dvDict

    def update(self, ptSetName, config=None):
        """This is the main routine for returning coordinates that have been
        updated by design variables. Multiple configs are not
        supported.

        Parameters
        ----------
        ptSetName : str
            Name of point-set to return. This must match ones of the
            given in an :func:`addPointSet()` call.
        """

        # get the new points
        newPts = numpy.zeros((self.points[ptSetName].nPts, 3))

        # we first need to update all points with their respective DVGeo objects
        for comp in self.compNames:
            ptsComp = self.comps[comp].DVGeo.update(ptSetName)

            # now save this info with the pointset mapping
            ptMap = self.points[ptSetName].compMap[comp]
            newPts[ptMap] = ptsComp

        # get the delta
        delta = newPts - self.points[ptSetName].points

        # then apply the intersection treatment
        for IC in self.intersectComps:
            # check if this IC is active for this ptSet
            if ptSetName in IC.points:
                delta = IC.update(ptSetName, delta)

        # now we are ready to take the delta which may be modified by the intersections
        newPts = self.points[ptSetName].points + delta

        # now, project the points that were warped back onto the trimesh
        for IC in self.intersectComps:
            if IC.projectFlag and ptSetName in IC.points:
                # new points will be modified in place using .... newPts array
                IC.project(ptSetName, newPts)

        # set the pointset up to date
        self.updated[ptSetName] = True

        return newPts

    def pointSetUpToDate(self, ptSetName):
        """
        This is used externally to query if the object needs to update
        its pointset or not. Essentially what happens, is when
        update() is called with a point set, it the self.updated dict
        entry for pointSet is flagged as true. Here we just return
        that flag. When design variables are set, we then reset all
        the flags to False since, when DVs are set, nothing (in
        general) will up to date anymore.

        Parameters
        ----------
        ptSetName : str
            The name of the pointset to check.
        """
        if ptSetName in self.updated:
            return self.updated[ptSetName]
        else:
            return True

    def getNDV(self):
        """Return the number of DVs"""
        # loop over components and sum number of DVs
        nDV = 0
        for comp in self.compNames:
            nDV += self.comps[comp].DVGeo.getNDV()
        return nDV

    def getVarNames(self):
        """
        Return a list of the design variable names. This is typically
        used when specifying a wrt= argument for pyOptSparse.

        Examples
        --------
        optProb.addCon(.....wrt=DVGeo.getVarNames())
        """
        dvNames = []
        # create a list of DVs from each comp
        for comp in self.compNames:
            # first get the list of DVs from this component
            varNames = self.comps[comp].DVGeo.getVarNames()

            for var in varNames:
                # then add the component's name to the DV name
                dvName = "%s:%s" % (comp, var)

                # finally append to the list
                dvNames.append(dvName)

        return dvNames

    def totalSensitivity(self, dIdpt, ptSetName, comm=None, config=None):
        """
        This function computes sensitivty information.

        Specificly, it computes the following:
        :math:`\\frac{dX_{pt}}{dX_{DV}}^T \\frac{dI}{d_{pt}}`

        Parameters
        ----------
        dIdpt : array of size (Npt, 3) or (N, Npt, 3)

            This is the total derivative of the objective or function
            of interest with respect to the coordinates in
            'ptSetName'. This can be a single array of size (Npt, 3)
            **or** a group of N vectors of size (Npt, 3, N). If you
            have many to do, it is faster to do many at once.

        ptSetName : str
            The name of set of points we are dealing with

        comm : MPI.IntraComm
            The communicator to use to reduce the final derivative. If
            comm is None, no reduction takes place.

        config : str or list
            Define what configurations this design variable will be applied to
            Use a string for a single configuration or a list for multiple
            configurations. The default value of None implies that the design
            variable appies to *ALL* configurations.


        Returns
        -------
        dIdxDict : dic
            The dictionary containing the derivatives, suitable for
            pyOptSparse

        Notes
        -----
        The ``child`` and ``nDVStore`` options are only used
        internally and should not be changed by the user.
        """
        # if comm:
        #     commPresent = True
        # else:
        #     commPresent = False

        # print('[%d] called totalSensitivity with comm:'%self.comm.rank, commPresent)

        # compute the total jacobian for this pointset
        # TODO, we dont even need to do this
        self._computeTotalJacobian(ptSetName)

        # compute IC jacobians if they are out of date
        # if not self.ICJupdated:
        #     self._computeICJacobian()

        # Make dIdpt at least 3D
        if len(dIdpt.shape) == 2:
            dIdpt = numpy.array([dIdpt])
        N = dIdpt.shape[0]

        # create a dictionary to save total sensitivity info that might come out of the ICs
        compSensList = []

        # if we projected points for any intersection treatment,
        # we need to propagate the derivative seed of the projected points
        # back to the seeds for the initial points we get after ID-warping
        for IC in self.intersectComps:
            if IC.projectFlag and ptSetName in IC.points:
                # initialize the seed contribution to the intersection seam and feature curves from project_b
                IC.seamBarProj[ptSetName] = numpy.zeros((N, IC.seam0.shape[0], IC.seam0.shape[1]))

                # we pass in dIdpt and the intersection object, along with pointset information
                # the intersection object adjusts the entries corresponding to projected points
                # and passes back dIdpt in place.
                compSens = IC.project_b(ptSetName, dIdpt, comm)

                # append this to the dictionary list...
                compSensList.append(compSens)

        # do the transpose multiplication

        # print('[%d] finished project_b'%self.comm.rank)

        # get the pointset
        # ptSet = self.points[ptSetName]

        # number of design variables
        # nDV = ptSet.jac.shape[1]

        # We should keep track of the intersections that this pointset is close to. There is no point in including the intersections far from this pointset in the sensitivity calc as the derivative seeds will be just zeros there.
        # ptSetICs = []

        # we need to go through all ICs bec even though some procs might not have points on the intersection,
        # communication is easier and we can reduce compSens as we compute them
        for IC in self.intersectComps:
            if ptSetName in IC.points:
                compSens = IC.sens(dIdpt, ptSetName, comm)
                # save the sensitivities from the intersection stuff
                compSensList.append(compSens)

            # instead this saves the dIdx for this intersection in the IC object
            # IC.sens(dIdpt, ptSetName, comm)

        # print('[%d] finished IC.sens'%self.comm.rank)

        # reshape the dIdpt array from [N] * [nPt] * [3] to  [N] * [nPt*3]
        dIdpt = dIdpt.reshape((dIdpt.shape[0], dIdpt.shape[1] * 3))

        # jacobian for the pointset
        jac = self.points[ptSetName].jac

        # this is the mat-vec product for the remaining seeds.
        # this only contains the effects of the FFD motion,
        # projections and intersections are handled separately in compSens
        dIdxT_local = jac.T.dot(dIdpt.T)
        dIdx_local = dIdxT_local.T

        # # now add the contributions from the intersections
        # for IC in self.intersectComps:
        #     dIdx_local = dIdx_local + IC.dIdx

        if comm:  # If we have a comm, globaly reduce with sum
            # print('[%d] before allreduce dIdx ='%self.comm.rank, dIdx_local)
            # comm.Barrier()
            dIdx = comm.allreduce(dIdx_local, op=MPI.SUM)
            # print('[%d] after  allreduce dIdx ='%self.comm.rank, dIdx_local)
            # comm.Barrier()
        else:
            dIdx = dIdx_local

        # print('[%d] finished comm.allreduce'%self.comm.rank)

        # use respective DVGeo's convert to dict functionality
        dIdxDict = OrderedDict()
        dvOffset = 0
        for comp in self.compNames:
            DVGeo = self.comps[comp].DVGeo
            nDVComp = DVGeo.getNDV()

            # print('[%d] full dIdx:'%(self.comm.rank), dIdx)
            # print('[%d] dIdx shape:'%self.comm.rank, dIdx.shape)

            # we only do this if this component has at least one DV
            if nDVComp > 0:
                # print('[%d] dIdx for comp %s:'%(self.comm.rank, comp), dIdx)
                # this part of the sensitivity matrix is owned by this dvgeo
                dIdxComp = DVGeo.convertSensitivityToDict(dIdx[:, dvOffset : dvOffset + nDVComp])

                # add the component names in front of the dictionary keys
                for k, v in dIdxComp.items():
                    dvName = "%s:%s" % (comp, k)
                    dIdxDict[dvName] = v

                # also increment the offset
                dvOffset += nDVComp

        # finally, we can add the contributions from triangulated component meshes
        for compSens in compSensList:
            # loop over the items of compSens, which are guaranteed to be in dIdxDict
            for k, v in compSens.items():
                # these will bring in effects from projections and intersection computations
                dIdxDict[k] += v

        # print('[%d] finished DVGeo.totalSensitivity!'%self.comm.rank)

        return dIdxDict

    def addVariablesPyOpt(
        self,
        optProb,
        globalVars=True,
        localVars=True,
        sectionlocalVars=True,
        ignoreVars=None,
        freezeVars=None,
        comps=None,
    ):
        """
        Add the current set of variables to the optProb object.

        Parameters
        ----------
        optProb : pyOpt_optimization class
            Optimization problem definition to which variables are added

        globalVars : bool
            Flag specifying whether global variables are to be added

        localVars : bool
            Flag specifying whether local variables are to be added

        ignoreVars : list of strings
            List of design variables the user DOESN'T want to use
            as optimization variables.

        freezeVars : list of string
            List of design variables the user WANTS to add as optimization
            variables, but to have the lower and upper bounds set at the current
            variable. This effectively eliminates the variable, but it the variable
            is still part of the optimization.

        comps : list of components we want to add DVs of. If no list is passed,
            we will add all dvs of all components
        """

        # if no list was provided, we use all components
        if comps is None:
            comps = self.compNames

        # we can simply loop over all DV objects and call their respective
        # addVariablesPyOpt function with the correct prefix.
        for comp in comps:
            self.comps[comp].DVGeo.addVariablesPyOpt(
                optProb,
                globalVars=globalVars,
                localVars=localVars,
                sectionlocalVars=sectionlocalVars,
                ignoreVars=ignoreVars,
                freezeVars=freezeVars,
                prefix=comp + ":",
            )

    def getLocalIndex(self, iVol, comp):
        """Return the local index mapping that points to the global
        coefficient list for a given volume"""
        # call this on the respective DVGeo
        DVGeo = self.comps[comp].DVGeo
        return DVGeo.FFD.topo.lIndex[iVol].copy()

    # ----------------------------------------------------------------------
    #        THE REMAINDER OF THE FUNCTIONS NEED NOT BE CALLED BY THE USER
    # ----------------------------------------------------------------------

    def _readCGNSFile(self, filename):
        # this function reads the unstructured CGNS grid in filename and returns
        # node coordinates and element connectivities.
        # Here, only the root proc reads the cgns file, broadcasts node and connectivity info.

        # create the featurecurve dictionary
        # curveConn = OrderedDict()

        # only root proc reads the file
        if self.comm.rank == 0:
            print("Reading file %s" % filename)
            # use the default routine in tsurftools
            nodes, sectionDict = tsurf_tools.getCGNSsections(filename, comm=MPI.COMM_SELF)
            print("Finished reading the cgns file")

            triConn = numpy.zeros((0, 3), dtype=numpy.int8)
            barsConn = {}
            # print('Part names in triangulated cgns file for %s'%filename)
            for part in sectionDict:
                # print(part)
                if "triaConnF" in sectionDict[part].keys():
                    # this is a surface, read the tri connectivities
                    triConn = numpy.vstack((triConn, sectionDict[part]["triaConnF"]))

                if "barsConn" in sectionDict[part].keys():
                    # this is a curve, save the curve connectivity
                    barsConn[part.lower()] = sectionDict[part]["barsConn"]

            print("The %s mesh has %d nodes and %d elements." % (filename, len(nodes), len(triConn)))
        else:
            # create these to recieve the data
            nodes = None
            triConn = None
            barsConn = None

        # each proc gets the nodes and connectivities
        # CHECK if this should be bcast or Bcast...
        nodes = self.comm.bcast(nodes, root=0)
        triConn = self.comm.bcast(triConn, root=0)
        barsConn = self.comm.bcast(barsConn, root=0)

        return nodes, triConn, barsConn

    def _computeTotalJacobian(self, ptSetName):
        """
        This routine computes the total jacobian. It takes the jacobians
        from respective DVGeo objects and also computes the jacobians for
        the intersection seams. We then use this information in the
        totalSensitivity function.
        """

        # number of design variables
        nDV = self.getNDV()

        # allocate space for the jacobian.
        jac = numpy.zeros((self.points[ptSetName].nPts * 3, nDV))

        # ptset
        ptSet = self.points[ptSetName]

        dvOffset = 0
        # we need to call computeTotalJacobian from all comps and get the jacobians for this pointset
        for comp in self.compNames:

            # number of design variables
            nDVComp = self.comps[comp].DVGeo.getNDV()

            # call the function to compute the total jacobian
            self.comps[comp].DVGeo.computeTotalJacobian(ptSetName)

            if self.comps[comp].DVGeo.JT[ptSetName] is not None:
                # we convert to dense storage.
                # this is probably not a good way to do this...
                # TODO: use sparse storage for these...
                compJ = self.comps[comp].DVGeo.JT[ptSetName].todense().T

                # loop over the entries and add one by one....
                # for i in range(len(ptSet.compMapFlat[comp])):
                #     jac[ptSet.compMapFlat[comp][i], dvOffset:dvOffset+nDVComp] = compJ[i,:]

                # or, do it (kinda) vectorized!
                jac[ptSet.compMapFlat[comp], dvOffset : dvOffset + nDVComp] = compJ[:, :]

            # increment the offset
            dvOffset += nDVComp

        # now we can save this jacobian in the pointset
        ptSet.jac = jac

    def _computeICJacobian(self):

        # loop over design variables and compute jacobians for all ICs

        # counts
        nDV = self.getNDV()
        nproc = self.comm.size
        rank = self.comm.rank

        # save the reference seams
        for IC in self.intersectComps:
            IC.seamRef = IC.seam.flatten()
            IC.jac = numpy.zeros((len(IC.seamRef), nDV))

        # We need to evaluate all the points on respective procs for FD computations

        # determine how many DVs this proc will perturb.
        n = 0
        for iDV in range(nDV):
            # I have to do this one.
            if iDV % nproc == rank:
                n += 1

        # perturb the DVs on different procs and compute the new point coordinates.
        # reqs = []
        for iDV in range(nDV):
            # I have to do this one.
            if iDV % nproc == rank:

                # TODO Use separate step sizes for different DVs?
                dh = self.dh

                # second order FD
                dvSave = self._getIthDV(iDV)
                self._setIthDV(iDV, dvSave + dh)

                # fwd step
                for IC in self.intersectComps:
                    IC.setSurface(MPI.COMM_SELF)
                    IC.seamPlus = IC.seam.copy()

                # bwd step and the actual FD comp
                self._setIthDV(iDV, dvSave - dh)
                for IC in self.intersectComps:
                    IC.setSurface(MPI.COMM_SELF)
                    IC.seamMinus = IC.seam.copy()

                    IC.jac[:, iDV] = (IC.seamPlus.flatten() - IC.seamMinus.flatten()) / (2 * dh)

                # first order FD
                # # Perturb the DV
                # dvSave = self._getIthDV(iDV)
                # self._setIthDV(iDV, dvSave+dh)

                # # Do any required intersections:
                # for IC in self.intersectComps:
                #     IC.setSurface(MPI.COMM_SELF)
                #     IC.jac[:,iDV] = (IC.seam.flatten() - IC.seamRef) / dh

                # Reset the DV
                self._setIthDV(iDV, dvSave)

        # Restore the seams
        for IC in self.intersectComps:
            IC.seam = IC.seamRef.reshape(len(IC.seam), 3)

        # also do the computation again on processors that perturbed the design
        for iDV in range(nDV):
            # I have to do this one.
            if iDV % nproc == rank:
                for IC in self.intersectComps:
                    # we only need to do this if we are using the projection for this intersection
                    # because we will need the intermediate values to be consistent when we are
                    # propagating the derivative seeds backward
                    if IC.projectFlag:
                        IC.setSurface(MPI.COMM_SELF)

        # loop over the DVs and scatter the perturbed points to original procs
        for iDV in range(nDV):

            # also bcast the intersection jacobians
            for IC in self.intersectComps:

                # create send/recv buffers
                if iDV % nproc == rank:
                    # we have to copy this because we need a contiguous array for bcast
                    buf = IC.jac[:, iDV].copy()
                else:
                    # receive buffer for procs that will get the result
                    buf = numpy.zeros(len(IC.seamRef))

                # bcast the intersection jacobian directly from the proc that perturbed this DV
                self.comm.Bcast([buf, len(IC.seamRef), MPI.DOUBLE], root=iDV % nproc)

                # set the value in the procs that dont have it
                if iDV % nproc != rank:
                    IC.jac[:, iDV] = buf.copy()

                # the slower version of the same bcast code
                # IC.jac[:,i] = self.comm.bcast(IC.jac[:,i], root=i%nproc)

        # set the flag
        self.ICJupdated = True

    def _setIthDV(self, iDV, val):
        # this function sets the design variable. The order is important, and we count every DV.

        nDVCum = 0
        # get the number of DVs from different DVGeos to figure out which DVGeo owns this DV
        for comp in self.compNames:
            # get the DVGeo object
            DVGeo = self.comps[comp].DVGeo
            # get the number of DVs on this DVGeo
            nDVComp = DVGeo.getNDV()
            # increment the cumulative number of DVs
            nDVCum += nDVComp

            # if we went past the index of the DV we are interested in
            if nDVCum > iDV:
                # this DVGeo owns the DV we want to set.
                xDV = DVGeo.getValues()

                # take back the global counter
                nDVCum -= nDVComp

                # loop over the DVs owned by this geo
                for k, v in xDV.items():
                    nDVKey = len(v)

                    # add the number of DVs for this DV key
                    nDVCum += nDVKey

                    # check if we went past the value we want to set
                    if nDVCum > iDV:
                        # take back the counter
                        nDVCum -= nDVKey

                        # this was an array...
                        if nDVKey > 1:
                            for i in range(nDVKey):
                                nDVCum += 1
                                if nDVCum > iDV:
                                    # finally...
                                    xDV[k][i] = val
                                    DVGeo.setDesignVars(xDV)
                                    return

                        # this is the value we want to set!
                        else:
                            xDV[k] = val
                            DVGeo.setDesignVars(xDV)
                            return

    def _getIthDV(self, iDV):
        # this function sets the design variable. The order is important, and we count every DV.

        nDVCum = 0
        # get the number of DVs from different DVGeos to figure out which DVGeo owns this DV
        for comp in self.compNames:
            # get the DVGeo object
            DVGeo = self.comps[comp].DVGeo
            # get the number of DVs on this DVGeo
            nDVComp = DVGeo.getNDV()
            # increment the cumulative number of DVs
            nDVCum += nDVComp

            # if we went past the index of the DV we are interested in
            if nDVCum > iDV:
                # this DVGeo owns the DV we want to set.
                xDV = DVGeo.getValues()

                # take back the global counter
                nDVCum -= nDVComp

                # loop over the DVs owned by this geo
                for k, v in xDV.items():
                    nDVKey = len(v)

                    # add the number of DVs for this DV key
                    nDVCum += nDVKey

                    # check if we went past the value we want to set
                    if nDVCum > iDV:
                        # take back the counter
                        nDVCum -= nDVKey

                        # this was an array...
                        if nDVKey > 1:
                            for i in range(nDVKey):
                                nDVCum += 1
                                if nDVCum > iDV:
                                    # finally...
                                    return xDV[k][i]

                        # this is the value we want to get!
                        else:
                            # assume single DVs come in arrays...
                            return xDV[k].copy()


class component:
    def __init__(self, name, DVGeo, nodes, triConn, barsConn, xMin, xMax):

        # save the info
        self.name = name
        self.DVGeo = DVGeo
        self.nodes = nodes
        self.triConn = triConn
        self.barsConn = barsConn
        self.xMin = xMin
        self.xMax = xMax

        # also a dictionary for DV names
        self.dvDict = {}

        # set a flag for triangulated meshes
        if nodes is None:
            self.triMesh = False
        else:
            self.triMesh = True

    def updateTriMesh(self):
        # update the triangulated surface mesh
        self.nodes = self.DVGeo.update("triMesh")


class PointSet:
    def __init__(self, points, comm):
        self.points = points
        self.nPts = len(self.points)
        self.compMap = OrderedDict()
        self.compMapFlat = OrderedDict()
        self.comm = comm


class CompIntersection:
    def __init__(
        self,
        compA,
        compB,
        dStarA,
        dStarB,
        featureCurves,
        distTol,
        DVGeo,
        project,
        marchDir,
        includeCurves,
        intDir,
        curveEpsDict,
        remeshBwd,
    ):
        """Class to store information required for an intersection.
        Here, we use some fortran code from pySurf.

        Input
        -----
        compA: ID of the first component
        compB: ID of the second component

        dStar, real : Radius over which to attenuate the deformation

        intDir, int : Direction of which intersection to pick,
                      +/- specifies the direction and
                      value (0,1,2) specifies the axis (x,y,z)

        featureCurves: list or dict. If user provides a list,
                       we pick the same marching direction for all
                       curves. if the user provides a dict, the keys
                       should be the curve names, and values should
                       be the march direction for each curve.

        Internally we will store the indices and the weights of the
        points that this intersection will have to modify. In general,
        all this code is not super efficient since it's all python
        """

        # same communicator with DVGeo
        self.comm = DVGeo.comm

        # counter for outputting curves etc at each update
        self.counter = 0

        # flag that determines if we will remesh the other side of the feature curves on compB
        self.remeshBwd = remeshBwd

        # tolerance used for each curve when mapping nodes to curves
        self.curveEpsDict = {}
        for k, v in curveEpsDict.items():
            self.curveEpsDict[k.lower()] = v

        # beginning and end element indices for each curve
        self.seamBeg = {}
        self.seamEnd = {}

        # indices of nodes to be projected to curves.
        self.curveProjIdx = {}

        # lists to track which feature curve is on which comp
        self.curvesOnA = []
        self.curvesOnB = []

        # a dict to to keep track which curves get points mapped to them
        # keys will be pointSetNames, then the values will be another dict,
        # where the keys are curve names and values will be a bool value
        self.curveProjFlag = {}

        # dicts to keep track of the coordinates and counts of points projected to curves
        self.curvePtCounts = {}
        self.curvePtCoords = {}

        # dict to keep track of the total number of points on each curve
        self.nCurvePts = {}

        # dictionary to keep the seam seeds that come from curve projections
        self.seamBarProj = {}

        # dictionaries to save the indices of points mapped to surfaces for each comp
        self.surfIdxA = {}
        self.surfIdxB = {}

        # names of compA and compB must be provided
        self.compA = DVGeo.comps[compA]
        self.compB = DVGeo.comps[compB]

        self.dStarA = dStarA
        self.dStarB = dStarB
        # self.halfdStar = self.dStar/2.0
        self.points = OrderedDict()

        # process the feature curves

        # list to save march directions
        marchDirs = []

        # list to save curve names where we remesh all the curve
        self.remeshAll = []

        # if a list is provided, we use this and the marchdir information
        if type(featureCurves) is list:
            self.featureCurveNames = featureCurves
            for i in range(len(self.featureCurveNames)):
                self.featureCurveNames[i] = self.featureCurveNames[i].lower()
                # get one march dir per curve
                marchDirs.append(marchDir)

        else:
            # if a dict is provided, the marchdirs are the dict values
            # we save this info in lists
            self.featureCurveNames = []
            # save the curve name and march direction information
            for k, v in featureCurves.items():
                self.featureCurveNames.append(k.lower())
                marchDirs.append(v)

        # now loop over the feature curves and flip if necessary
        for ii, curveName in enumerate(self.featureCurveNames):
            # figure out which comp owns this curve...
            if curveName in self.compB.barsConn:
                curveComp = self.compB
                self.curvesOnB.append(curveName)
            elif curveName in self.compA.barsConn:
                curveComp = self.compA
                self.curvesOnA.append(curveName)
            else:
                raise Error("Curve %s does not belong in %s or %s" % (curveName, self.compA.name, self.compB.name))

            # sort the feature curve
            newConn, newMap = tsurf_tools.FEsort(curveComp.barsConn[curveName].tolist())

            # we only want to have a single curve
            if len(newConn) > 1:
                raise Error("the curve %s generated more than one curve with FESort" % curveName)

            # get the connectivity
            newConn = newConn[0]

            # we may also need to flip the curve
            curveNodes = curveComp.nodes

            if marchDirs[ii] is None:
                # we remesh all of this curve
                self.remeshAll.append(curveName)
                curveComp.barsConn[curveName] = newConn
            else:
                # get the direction we want to march
                mdir = abs(marchDirs[ii]) - 1
                msign = numpy.sign(marchDirs[ii])

                # check if we need to flip
                if msign * curveNodes[newConn[0][0]][mdir] > msign * curveNodes[newConn[0][1]][mdir]:
                    # flip on both axes bec. pleiades complains when we flip both at the same time
                    newConn = numpy.flip(newConn, axis=0)
                    newConn = numpy.flip(newConn, axis=1)

                # save the new connectivity
                curveComp.barsConn[curveName] = newConn

        self.distTol = distTol

        # flag to determine if we want to project nodes after intersection treatment
        self.projectFlag = project

        # create the dictionary if we are projecting.
        if project:
            self.projData = {}
            if includeCurves:
                # dict to save the data related to projection to curves
                self.curveProjData = {}

        # flag to include feature curves in ID-warping
        self.incCurves = includeCurves

        # direction to pick if we have multiple intersection curves
        self.intDir = intDir

        # only the node coordinates will be modified for the intersection calculations because we have calculated and saved all the connectivity information
        if self.comm.rank == 0:
            print("Computing initial intersection between %s and %s" % (compA, compB))
        self.seam0 = self._getIntersectionSeam(self.comm, firstCall=True)
        self.seam = self.seam0.copy()

    def setSurface(self, comm):
        """This set the new udpated surface on which we need to compute the new intersection curve"""

        # get the updated surface coordinates
        self._getUpdatedCoords()

        self.seam = self._getIntersectionSeam(comm)

    def addPointSet(self, pts, ptSetName, compMap, comm):

        # Figure out which points this intersection object has to deal with.

        # use Ney's fortran code to project the point on curve
        # Get number of points
        nPoints = len(pts)

        # Initialize references if user provided none
        dist2 = numpy.ones(nPoints) * 1e10
        xyzProj = numpy.zeros((nPoints, 3))
        tanProj = numpy.zeros((nPoints, 3))
        elemIDs = numpy.zeros((nPoints), dtype="int32")

        # only call the fortran code if we have at least one point
        if nPoints > 0:
            # Call fortran code
            # This will modify xyzProj, tanProj, dist2, and elemIDs if we find better projections than dist2.
            # Remember that we should adjust some indices before calling the Fortran code
            # Remember to use [:] to don't lose the pointer (elemIDs is an input/output variable)
            elemIDs[:] = (
                elemIDs + 1
            )  # (we need to do this separetely because Fortran will actively change elemIDs contents.
            curveSearchAPI.curvesearchapi.mindistancecurve(
                pts.T, self.nodes0.T, self.conn0.T + 1, xyzProj.T, tanProj.T, dist2, elemIDs
            )

            # Adjust indices back to Python standards
            elemIDs[:] = elemIDs - 1

        # dist2 has the array of squared distances
        d = numpy.sqrt(dist2)

        indices = []
        factors = []
        for i in range(len(pts)):
            # figure out which component this point is mapped to
            if i in compMap[self.compA.name]:
                # component A owns this
                dStar = self.dStarA
            else:
                # comp B owns this point
                dStar = self.dStarB

            # then get the halfdStar for that component
            halfdStar = dStar / 2.0

            if d[i] < dStar:

                # Compute the factor
                if d[i] < halfdStar:
                    factor = 0.5 * (d[i] / halfdStar) ** 3
                else:
                    factor = 0.5 * (2 - ((dStar - d[i]) / halfdStar) ** 3)

                # Save the index and factor
                indices.append(i)
                factors.append(factor)

        # # Print the number of points that get associated with the intersection.
        # nPointGlobal = self.comm.allreduce(len(factors), op=MPI.SUM)
        # if self.comm.rank == 0:
        #     intName = vsp.GetContainerName(self.compA)+'_'+vsp.GetContainerName(self.compB)
        #     print('DVGEO VSP:\n%d points associated with intersection %s'%(nPointGlobal, intName))

        # Save the affected indices and the factor in the little dictionary
        self.points[ptSetName] = [pts.copy(), indices, factors, comm]

        # now we need to figure out which components we are projecting to if projection is enabled
        # this can be done faster above but whatever
        if self.projectFlag:
            flagA = False
            flagB = False

            indices = self.points[ptSetName][1]

            # create the list we use to map the points to projection components
            indA = []
            indB = []

            # maybe we can do this vectorized
            for ind in indices:

                # check compA
                if ind in compMap[self.compA.name]:
                    flagA = True
                    indA.append(ind)

                # check compB
                if ind in compMap[self.compB.name]:
                    flagB = True
                    indB.append(ind)

            # now we create the dictionaries to propogate projection data
            self.projData[ptSetName] = {
                # we need two dictionaries.
                # we will save more data here, but for now, just save the flag and indices
                "compA": {"flag": flagA, "ind": indA},
                "compB": {"flag": flagB, "ind": indB},
            }

            # if we include the feature curves in the warping, we also need to project the added points to the intersection and feature curves and determine how the points map to the curves
            if self.incCurves:

                # convert the list to an array
                # we specify the dtype because numpy cannot know the type when 'indices' is empty
                indices = numpy.array(indices, dtype="intc")

                # get the coordinates of all points affected by this intersection
                ptsToCurves = pts[indices]

                # project these to the combined curves
                # use Ney's fortran code to project the point on curve
                # Get number of points
                nPoints = len(ptsToCurves)

                # Initialize references if user provided none
                dist2 = numpy.ones(nPoints) * 1e10
                xyzProj = numpy.zeros((nPoints, 3))
                tanProj = numpy.zeros((nPoints, 3))
                elemIDs = numpy.zeros((nPoints), dtype="int32")

                # only call the fortran code if we have at least one point
                if nPoints > 0:
                    # Call fortran code
                    # This will modify xyzProj, tanProj, dist2, and elemIDs if we find better projections than dist2.
                    # Remember that we should adjust some indices before calling the Fortran code
                    # Remember to use [:] to don't lose the pointer (elemIDs is an input/output variable)
                    elemIDs[:] = elemIDs + 1
                    # (we need to do this separetely because Fortran will actively change elemIDs contents.
                    curveSearchAPI.curvesearchapi.mindistancecurve(
                        ptsToCurves.T, self.seam0.T, self.seamConn.T + 1, xyzProj.T, tanProj.T, dist2, elemIDs
                    )

                    # Adjust indices back to Python standards
                    elemIDs[:] = elemIDs - 1

                # dist2 has the array of squared distances
                d = numpy.sqrt(dist2)

                # get the names of all curves including the intersection
                allCurves = ["intersection"]
                for curveName in self.featureCurveNames:
                    allCurves.append(curveName)

                # track the points that dont get associated with any curve
                # get a full masking array with zeros
                allNodesBool = numpy.zeros(len(elemIDs))

                # dict to save the pt indices
                self.curveProjIdx[ptSetName] = {}
                # dict to save other data
                self.curveProjData[ptSetName] = {}

                # now loop over feature curves and use the epsilon that each curve has
                # to determine which points maps to which curve
                for curveName in allCurves:

                    # get the epsilon for this curve
                    # we map the points closer than eps to this curve
                    eps = self.curveEpsDict[curveName]

                    # also get the range of element IDs this curve owns
                    seamBeg = self.seamBeg[curveName]
                    seamEnd = self.seamEnd[curveName]

                    # this returns a bool array of indices that satisfy the conditions
                    # we check for elemIDs because we projected to all curves at once
                    curveBool = numpy.all([d < eps, elemIDs >= seamBeg, elemIDs < seamEnd], axis=0)

                    # get the indices of the points mapped to this element
                    idxs = numpy.nonzero(curveBool)

                    # save the indices. idx has the indices of the "indices" array
                    # that had the indices of points that get any intersection treatment
                    # print("[%d] curvename: %s indices[idxs] idxs"%(self.comm.rank, curveName), indices[idxs], idxs)
                    self.curveProjIdx[ptSetName][curveName] = numpy.array(indices[idxs])

                    # uncomment to get the output
                    # ptCoords = ptsToCurves[idxs]
                    # tecplot_interface.write_tecplot_scatter('%s.plt'%curveName, curveName, ['X', 'Y', 'Z'], ptCoords)

                    # also update the masking array
                    # we will use this to figure out the indices that did not get attached to any curves
                    allNodesBool = numpy.any([curveBool, allNodesBool], axis=0)

                    # create an empty dict to save aux data later on
                    self.curveProjData[ptSetName][curveName] = {}

                # negate the surface mask and get indices
                surfPtIdx = numpy.nonzero(numpy.logical_not(allNodesBool))

                # figure out which of these surfNodes live only on components A and B
                allSurfIdx = numpy.array(indices[surfPtIdx[0]])

                # component A
                mask = numpy.in1d(allSurfIdx, indA, assume_unique=True)
                # this is the local indices of the points affected
                self.surfIdxA[ptSetName] = allSurfIdx[numpy.nonzero(mask)]

                # component B
                mask = numpy.in1d(allSurfIdx, indB, assume_unique=True)
                self.surfIdxB[ptSetName] = allSurfIdx[numpy.nonzero(mask)]

                # initialize the bool dict for this pointset
                self.curveProjFlag[ptSetName] = {}
                self.curvePtCounts[ptSetName] = {}
                self.curvePtCoords[ptSetName] = {}
                self.nCurvePts[ptSetName] = {}

                # we need to figure out if we have any points mapped to curves on comp A
                for curveName in allCurves:

                    # initialize the flag to false
                    # self.curveProjFlag[ptSetName][curveName] = False

                    # get the indices mapped to this curve, on this proc
                    idxs = self.curveProjIdx[ptSetName][curveName]

                    # call the utility function
                    nPtsTotal, nPtsProcs, curvePtCoords = self._commCurveProj(pts, idxs, comm)

                    # for i in range(comm.size):
                    # if i == comm.rank:
                    # print('[%d] data for %s:'%(self.comm.rank, curveName), '\nnptsTotal:',nPtsTotal, '\nnPtsProcs:', nPtsProcs,'\n',flush=True)
                    # comm.Barrier()

                    # check if we have at least one point globally that mapped to this curve
                    # if nPtsTotal > 0:
                    # set the flag, we do have projected points on this curve
                    # self.curveProjFlag[ptSetName][curveName] = True

                    # save the displacements and points
                    self.curvePtCounts[ptSetName][curveName] = nPtsProcs
                    self.curvePtCoords[ptSetName][curveName] = curvePtCoords

                    # also save the total number for convenience
                    self.nCurvePts[ptSetName][curveName] = nPtsTotal

    def update(self, ptSetName, delta):

        """Update the delta in ptSetName with our correction. The delta need
        to be supplied as we will be changing it and returning them
        """

        # original coordinates of the added pointset
        pts = self.points[ptSetName][0]
        # indices of the points that get affected by this intersection
        indices = self.points[ptSetName][1]
        # factors for each node in pointSet
        factors = self.points[ptSetName][2]

        # get the comm
        # comm = self.points[ptSetName][3]
        # if comm:
        #     rank = comm.rank
        # else:
        #     rank = 0

        # coordinates for the remeshed curves
        # we use the initial seam coordinates here
        coor = self.seam0
        # bar connectivity for the remeshed elements
        conn = self.seamConn
        # deltas for each point (nNode, 3) in size
        dr = self.seam - self.seam0

        # define an epsilon to avoid dividing by zero later on
        eps = 1e-50  # 1e-32

        # loop over the points that get affected
        for i in range(len(factors)):
            # j is the index of the point in the full set we are working with.
            j = indices[i]

            # coordinates of the original point
            rp = pts[j]

            # Run the weighted interp:
            # num = numpy.zeros(3)
            # den = 0.0

            # # we loop over elements and compute the elemnent-wise integral on each line
            # for k in range(len(conn)):
            #     # get the two end points for this line
            #     r0 = coor[conn[k,0]]
            #     r1 = coor[conn[k,1]]

            #     # get the deltas for two end points
            #     dr0 = dr[conn[k,0]]
            #     dr1 = dr[conn[k,1]]

            #     # compute a, b, and c coefficients
            #     a = (r1[0]-r0[0])**2 + (r1[1]-r0[1])**2 + (r1[2]-r0[2])**2
            #     b = 2 * ((r1[0]-r0[0])*(r0[0]-rp[0]) + (r1[1]-r0[1])*(r0[1]-rp[1]) + (r1[2]-r0[2])*(r0[2]-rp[2]))
            #     c = (r0[0]-rp[0])**2 + (r0[1]-rp[1])**2 + (r0[2]-rp[2])**2

            #     # compute some re-occurring terms
            #     det = b*b - 4*a*c + 1e-32 # add an epsilon so that the determinant never becomes zero
            #     sabc = numpy.sqrt(a+b+c)
            #     sc   = numpy.sqrt(c)

            #     # denominators on the integral evaluations
            #     den1 = det*sabc
            #     den2 = det*sc

            #     # integral evaluations
            #     eval1 = -2*(2*a +b)/den1 + 2*b/den2
            #     eval2 = (2*b + 4*c)/den1 - 4*c/den2

            #     # numerator gets two integrals with the delta components
            #     num += (dr1-dr0) * eval2 + dr0 * eval1

            #     # denominator only gets one integral
            #     den += eval1

            # interp = num / den

            # Do it vectorized!

            # get the two end points for the line elements
            r0 = coor[conn[:, 0]]
            r1 = coor[conn[:, 1]]

            # get the deltas for two end points
            dr0 = dr[conn[:, 0]]
            dr1 = dr[conn[:, 1]]

            # compute a, b, and c coefficients
            a = (r1[:, 0] - r0[:, 0]) ** 2 + (r1[:, 1] - r0[:, 1]) ** 2 + (r1[:, 2] - r0[:, 2]) ** 2
            b = 2 * (
                (r1[:, 0] - r0[:, 0]) * (r0[:, 0] - rp[0])
                + (r1[:, 1] - r0[:, 1]) * (r0[:, 1] - rp[1])
                + (r1[:, 2] - r0[:, 2]) * (r0[:, 2] - rp[2])
            )
            c = (r0[:, 0] - rp[0]) ** 2 + (r0[:, 1] - rp[1]) ** 2 + (r0[:, 2] - rp[2]) ** 2

            # distances for each element
            dists = numpy.sqrt(numpy.maximum(a, 0.0))

            # compute some re-occurring terms
            # the determinant can be zero or negative, but it CANNOT be positive
            # this is because the quadratic that defines the distance from the line cannot have two roots.
            # if the point is on the line, the quadratic will have a single root...
            det = b * b - 4 * a * c
            # these might be negative 1e-20sth so clip them...
            # these will be strictly zero or greater than zero.
            # Numerically, they cannot be negative bec. we are working with real numbers
            sabc = numpy.sqrt(numpy.maximum(a + b + c, 0.0))
            sc = numpy.sqrt(numpy.maximum(c, 0.0))

            # denominators on the integral evaluations
            # add an epsilon so that these terms never become zero
            # det <= 0, sabc and sc >= 0, therefore the den1 and den2 should be <=0
            den1 = det * sabc - eps
            den2 = det * sc - eps

            # integral evaluations
            eval1 = (-2 * (2 * a + b) / den1 + 2 * b / den2) * dists
            eval2 = ((2 * b + 4 * c) / den1 - 4 * c / den2) * dists

            # denominator only gets one integral
            den = numpy.sum(eval1)

            # do each direction separately
            interp = numpy.zeros(3)
            for iDim in range(3):
                # numerator gets two integrals with the delta components
                num = numpy.sum((dr1[:, iDim] - dr0[:, iDim]) * eval2 + dr0[:, iDim] * eval1)
                # final result
                interp[iDim] = num / den

            # Now the delta is replaced by 1-factor times the weighted
            # interp of the seam * factor of the original:
            delta[j] = factors[i] * delta[j] + (1 - factors[i]) * interp
            # delta[j] = interp

        # if comm:
        #     comm.Barrier()

        return delta

    def update_d(self, ptSetName, dPt, dSeam):

        """forward mode differentiated version of the update routine.
        Note that dPt and dSeam are both one dimensional arrays
        """
        pts = self.points[ptSetName][0]
        indices = self.points[ptSetName][1]
        factors = self.points[ptSetName][2]

        # we need to reshape the arrays for simpler code
        dSeam = dSeam.reshape((len(self.seam0), 3))

        for i in range(len(factors)):
            # j is the index of the point in the full set we are
            # working with.
            j = indices[i]

            # Do it vectorized
            rr = pts[j] - self.seam0
            LdefoDist = 1.0 / numpy.sqrt(rr[:, 0] ** 2 + rr[:, 1] ** 2 + rr[:, 2] ** 2 + 1e-16)
            LdefoDist3 = LdefoDist ** 3
            Wi = LdefoDist3
            den = numpy.sum(Wi)
            interp_d = numpy.zeros(3)
            for iDim in range(3):
                interp_d[iDim] = numpy.sum(Wi * dSeam[:, iDim]) / den

                # Now the delta is replaced by 1-factor times the weighted
                # interp of the seam * factor of the original:
                dPt[j * 3 + iDim] = factors[i] * dPt[j * 3 + iDim] + (1 - factors[i]) * interp_d[iDim]

        return

    def sens(self, dIdPt, ptSetName, comm):
        # Return the reverse accumulation of dIdpt on the seam
        # nodes. Also modifies the dIdp array accordingly.

        # original coordinates of the added pointset
        pts = self.points[ptSetName][0]
        # indices of the points that get affected by this intersection
        indices = self.points[ptSetName][1]
        # factors for each node in pointSet
        factors = self.points[ptSetName][2]

        # coordinates for the remeshed curves
        # we use the initial seam coordinates here
        coor = self.seam0
        # bar connectivity for the remeshed elements
        conn = self.seamConn

        # define an epsilon to avoid dividing by zero later on
        eps = 1e-50

        # if we are handling more than one function,
        # seamBar will contain the seeds for each function separately
        seamBar = numpy.zeros((dIdPt.shape[0], self.seam0.shape[0], self.seam0.shape[1]))

        # if we have the projection flag, then we need to add the contribution to seamBar from that
        if self.projectFlag:
            seamBar += self.seamBarProj[ptSetName]

        # TODO we can vectorize the k-loop
        for k in range(dIdPt.shape[0]):
            for i in range(len(factors)):

                # j is the index of the point in the full set we are working with.
                j = indices[i]

                # coordinates of the original point
                rp = pts[j]

                # This is the local seed (well the 3 seeds for the point)
                localVal = dIdPt[k, j, :] * (1 - factors[i])

                # Scale the dIdpt by the factor..dIdpt is input/output
                dIdPt[k, j, :] *= factors[i]

                # get the two end points for the line elements
                r0 = coor[conn[:, 0]]
                r1 = coor[conn[:, 1]]
                # compute a, b, and c coefficients
                a = (r1[:, 0] - r0[:, 0]) ** 2 + (r1[:, 1] - r0[:, 1]) ** 2 + (r1[:, 2] - r0[:, 2]) ** 2
                b = 2 * (
                    (r1[:, 0] - r0[:, 0]) * (r0[:, 0] - rp[0])
                    + (r1[:, 1] - r0[:, 1]) * (r0[:, 1] - rp[1])
                    + (r1[:, 2] - r0[:, 2]) * (r0[:, 2] - rp[2])
                )
                c = (r0[:, 0] - rp[0]) ** 2 + (r0[:, 1] - rp[1]) ** 2 + (r0[:, 2] - rp[2]) ** 2
                # distances for each element
                dists = numpy.sqrt(numpy.maximum(a, 0.0))

                # compute some re-occurring terms
                det = b * b - 4 * a * c
                sabc = numpy.sqrt(numpy.maximum(a + b + c, 0.0))
                sc = numpy.sqrt(numpy.maximum(c, 0.0))
                # denominators on the integral evaluations
                den1 = det * sabc - eps
                den2 = det * sc - eps
                # integral evaluations
                eval1 = (-2 * (2 * a + b) / den1 + 2 * b / den2) * dists
                eval2 = ((2 * b + 4 * c) / den1 - 4 * c / den2) * dists

                # denominator only gets one integral
                den = numpy.sum(eval1)

                evalDiff = eval1 - eval2

                # do each direction separately
                for iDim in range(3):
                    # seeds for the r0 point
                    seamBar[k, conn[:, 0], iDim] += localVal[iDim] * evalDiff / den

                    # seeds for the r1 point
                    seamBar[k, conn[:, 1], iDim] += localVal[iDim] * eval2 / den

        # seamBar is the bwd seeds for the intersection curve...
        # it is N,nseampt,3 in size
        # print('[%d] calling getIntersectionSeam_b'%comm.rank)
        # now call the reverse differentiated seam computation
        compSens = self._getIntersectionSeam_b(seamBar, comm)
        # print('[%d] after getIntersectionSeam_b'%comm.rank)

        # instead, multiply this with the transpose of the jacobian
        # seamBar = seamBar.reshape( [dIdPt.shape[0], 3*self.seam0.shape[0]] )

        # sb = numpy.zeros((dIdPt.shape[0],3*self.seam0.shape[0] ))
        # for i in range(dIdPt.shape[0]):
        #     sb[i,:] = seamBar[i,:,:].flatten()

        # dIdx = ((self.jac.T).dot(sb.T)).T

        # print(dIdx.shape)
        # quit()

        # now convert dIdx to dictionary

        # self.dIdx = dIdx

        return compSens

    def project(self, ptSetName, newPts):
        # we need to build ADTs for both components if we have any components that lie on either
        # we also need to save ALL intermediate variables for gradient computations in reverse mode

        # get the comm for this point set
        comm = self.points[ptSetName][3]

        # if comm:
        #     rank = comm.rank
        # else:
        #     rank = 0

        self.comm.Barrier()

        # check if we need to worry about either surface
        # we will use these flags to figure out if we need to do warping.
        # we need to do the comm for the updated curves regardless
        flagA = False
        flagB = False
        if len(self.curvesOnA) > 0:
            # if len(self.surfIdxA[ptSetName]) > 0 and len(self.curvesOnA) > 0:
            flagA = True
        if len(self.curvesOnB) > 0:
            # if len(self.surfIdxB[ptSetName]) > 0 and len(self.curvesOnB) > 0:
            flagB = True

        # do the pts on the intersection outside the loop
        nptsg = self.nCurvePts[ptSetName]["intersection"]

        # the deltas for these points are zero. they should already be on the intersection
        # also get the initial coordinates of these points. We use this during warping
        # intersection curve will be on both components
        if flagA:
            deltaA = numpy.zeros((nptsg, 3))
            curvePtCoordsA = self.curvePtCoords[ptSetName]["intersection"].copy()
        else:
            deltaA = numpy.zeros((0, 3))
            curvePtCoordsA = numpy.zeros((0, 3))

        if flagB:
            deltaB = numpy.zeros((nptsg, 3))
            curvePtCoordsB = self.curvePtCoords[ptSetName]["intersection"].copy()
        else:
            deltaB = numpy.zeros((0, 3))
            curvePtCoordsB = numpy.zeros((0, 3))

        # tecplot_interface.write_tecplot_scatter('intersection_warped_pts.plt', 'intersection', ['X', 'Y', 'Z'], newPts[idx])

        # loop over the feature curves that we need to project
        for curveName in self.featureCurveNames:

            # get the indices of points we need to project
            idx = self.curveProjIdx[ptSetName][curveName]
            # print('[%d] curvename: %s idx'%(rank, curveName), idx.shape, idx)

            # these are the updated coordinates that will be projected to the curve
            ptsOnCurve = newPts[idx, :].copy()

            # tecplot_interface.write_tecplot_scatter('%s_warped_pts.plt'%curveName, curveName, ['X', 'Y', 'Z'], ptsOnCurve)

            # conn of the current curve
            seamBeg = self.seamBeg[curveName]
            seamEnd = self.seamEnd[curveName]
            curveConn = self.seamConn[seamBeg:seamEnd]

            # project these to the combined curves
            # use Ney's fortran code to project the point on curve
            # Get number of points
            nPoints = ptsOnCurve.shape[0]
            # print('[%d] curveName: %s, nPoints on the fwd pass: %d'%(self.comm.rank, curveName, nPoints))

            # Initialize references if user provided none
            dist2 = numpy.ones(nPoints) * 1e10
            xyzProj = numpy.zeros((nPoints, 3))
            tanProj = numpy.zeros((nPoints, 3))
            elemIDs = numpy.zeros((nPoints), dtype="int32")

            # only call the fortran code if we have at least one point
            if nPoints > 0:
                # Call fortran code
                # This will modify xyzProj, tanProj, dist2, and elemIDs if we find better projections than dist2.
                # Remember that we should adjust some indices before calling the Fortran code
                # Remember to use [:] to don't lose the pointer (elemIDs is an input/output variable)
                elemIDs[:] = (
                    elemIDs + 1
                )  # (we need to do this separetely because Fortran will actively change elemIDs contents.
                curveMask = curveSearchAPI.curvesearchapi.mindistancecurve(
                    ptsOnCurve.T, self.seam.T, curveConn.T + 1, xyzProj.T, tanProj.T, dist2, elemIDs
                )

                # Adjust indices back to Python standards
                elemIDs[:] = elemIDs - 1

                # we only have the curvemask if we do the projection on this proc
                self.curveProjData[ptSetName][curveName]["curveMask"] = curveMask

            # dist2 has the array of squared distances
            # d = numpy.sqrt(dist2)

            # save some information for gradient comp
            self.curveProjData[ptSetName][curveName]["xyz"] = ptsOnCurve.copy()
            self.curveProjData[ptSetName][curveName]["coor"] = self.seam.copy()
            self.curveProjData[ptSetName][curveName]["barsConn"] = curveConn.copy()
            self.curveProjData[ptSetName][curveName]["xyzProj"] = xyzProj.copy()
            self.curveProjData[ptSetName][curveName]["tanProj"] = tanProj.copy()
            self.curveProjData[ptSetName][curveName]["elemIDs"] = elemIDs.copy()

            # get the delta for the points on this proc
            deltaLocal = xyzProj - ptsOnCurve

            # tecplot_interface.write_tecplot_scatter('%s_projected_pts.plt'%curveName, curveName, ['X', 'Y', 'Z'], xyzProj)

            # update the point coordinates on this processor.
            # we do not need to do any communication for this
            # since newPts is the final coordinates of the points we just projected
            newPts[idx] = xyzProj

            # communicate the deltas
            if comm:
                sizes = self.curvePtCounts[ptSetName][curveName]
                disp = numpy.array([numpy.sum(sizes[:i]) for i in range(comm.size)], dtype="intc")

                # save these for grad comp
                self.curveProjData[ptSetName][curveName]["sizes"] = sizes
                self.curveProjData[ptSetName][curveName]["disp"] = disp

                # sendbuf
                deltaLocal = deltaLocal.flatten()
                sendbuf = [deltaLocal, sizes[comm.rank] * 3]

                # recvbuf
                nptsg = self.nCurvePts[ptSetName][curveName]
                deltaGlobal = numpy.zeros(nptsg * 3)
                recvbuf = [deltaGlobal, sizes * 3, disp * 3, MPI.DOUBLE]

                # do an allgatherv
                comm.Allgatherv(sendbuf, recvbuf)

                # reshape into a nptsg,3 array
                deltaGlobal = deltaGlobal.reshape((nptsg, 3))

            else:
                # we dont have a comm, so this is a "serial" pointset
                deltaGlobal = deltaLocal

                # also save the sizes and disp stuff as if we have one proc
                self.curveProjData[ptSetName][curveName]["sizes"] = self.curvePtCounts[ptSetName][curveName]
                self.curveProjData[ptSetName][curveName]["disp"] = [0]

            # we only add the deltaLocal to deltaA if this curve is on compA,
            # and we have points on compA surface
            if curveName in self.curvesOnA and flagA:
                # stack the deltas
                deltaA = numpy.vstack((deltaA, deltaGlobal))

                # also stack the original coordinates for warping
                curvePtCoordsNew = self.curvePtCoords[ptSetName][curveName]
                curvePtCoordsA = numpy.vstack((curvePtCoordsA, curvePtCoordsNew))

            # do the same for compB
            # we can use elif bec. one curve cannot be on both comps
            elif curveName in self.curvesOnB and flagB:
                # stack the deltas
                deltaB = numpy.vstack((deltaB, deltaGlobal))

                # also stack the original coordinates for warping
                curvePtCoordsNew = self.curvePtCoords[ptSetName][curveName]
                curvePtCoordsB = numpy.vstack((curvePtCoordsB, curvePtCoordsNew))

        self.comm.Barrier()

        # then, we warp all of the nodes that were affected by the intersection treatment, using the deltas from the previous project to curve step

        if flagA:
            self._warpSurfPts(self.points[ptSetName][0], newPts, self.surfIdxA[ptSetName], curvePtCoordsA, deltaA)

        if flagB:
            self._warpSurfPts(self.points[ptSetName][0], newPts, self.surfIdxB[ptSetName], curvePtCoordsB, deltaB)

        # save some info for the sens. computations
        self.curveProjData[ptSetName]["curvePtCoordsA"] = curvePtCoordsA
        self.curveProjData[ptSetName]["curvePtCoordsB"] = curvePtCoordsB

        # if comm:
        # comm.Barrier()

        # get the flags for components
        flagA = self.projData[ptSetName]["compA"]["flag"]
        flagB = self.projData[ptSetName]["compB"]["flag"]

        # call the actual driver with the info to prevent code multiplication
        if flagA:
            # get the indices of points that get projected to compA
            indA = self.projData[ptSetName]["compA"]["ind"]
            # get the points using the mapping
            ptsA = newPts[indA]
            # call the projection routine with the info
            # this returns the projected points and we use the same mapping to put them back in place
            newPts[indA] = self._projectToComponent(ptsA, self.compA, self.projData[ptSetName]["compA"])

        # do the same for B
        if flagB:
            indB = self.projData[ptSetName]["compB"]["ind"]
            ptsB = newPts[indB]
            newPts[indB] = self._projectToComponent(ptsB, self.compB, self.projData[ptSetName]["compB"])

    def project_b(self, ptSetName, dIdpt, comm):
        # call the functions to propagate ad seeds bwd
        # we need to build ADTs for both components if we have any components that lie on either
        # we also need to save ALL intermediate variables for gradient computations in reverse mode

        # number of functions we have
        N = dIdpt.shape[0]

        # get the flags for components
        flagA = self.projData[ptSetName]["compA"]["flag"]
        flagB = self.projData[ptSetName]["compB"]["flag"]

        # dictionary to accumulate triangulated mesh sensitivities
        compSens_local = {}

        # call the actual driver with the info to prevent code multiplication
        if flagA:
            # get the indices of points that get projected to compA
            indA = self.projData[ptSetName]["compA"]["ind"]
            # get the points using the mapping
            dIdptA = dIdpt[:, indA]
            # call the projection routine with the info
            # this returns the projected points and we use the same mapping to put them back in place
            dIdpt[:, indA], compSensA = self._projectToComponent_b(
                dIdptA, self.compA, self.projData[ptSetName]["compA"]
            )

            for k, v in compSensA.items():
                kNew = "%s:%s" % (self.compA.name, k)
                compSens_local[kNew] = v

        # set the compSens entries to all zeros on these procs
        else:
            # get the values from each DVGeo
            xA = self.compA.DVGeo.getValues()

            # loop over each entry in xA and xB and create a dummy zero gradient array for all
            for k, v in xA.items():
                kNew = "%s:%s" % (self.compA.name, k)
                # create the zero array:
                zeroSens = numpy.zeros((N, v.shape[0]))
                compSens_local[kNew] = zeroSens

        # do the same for B
        if flagB:
            indB = self.projData[ptSetName]["compB"]["ind"]
            dIdptB = dIdpt[:, indB]
            dIdpt[:, indB], compSensB = self._projectToComponent_b(
                dIdptB, self.compB, self.projData[ptSetName]["compB"]
            )

            for k, v in compSensB.items():
                kNew = "%s:%s" % (self.compB.name, k)
                compSens_local[kNew] = v
        # set the compSens entries to all zeros on these procs
        else:
            # get the values from each DVGeo
            xB = self.compB.DVGeo.getValues()

            # loop over each entry in xA and xB and create a dummy zero gradient array for all
            for k, v in xB.items():
                kNew = "%s:%s" % (self.compB.name, k)
                # create the zero array:
                zeroSens = numpy.zeros((N, v.shape[0]))
                compSens_local[kNew] = zeroSens

        # finally sum the results across procs if we are provided with a comm
        if comm:
            compSens = {}
            # because the results are in a dictionary, we need to loop over the items and sum
            for k in compSens_local:
                compSens[k] = comm.allreduce(compSens_local[k], op=MPI.SUM)
        else:
            # we can just pass the dictionary
            compSens = compSens_local

        # now we do the warping

        # get the comm for this point set
        ptSetComm = self.points[ptSetName][3]

        # check if any processor did any warping on compA
        curvePtCoordsA = self.curveProjData[ptSetName]["curvePtCoordsA"]
        curvePtCoordsB = self.curveProjData[ptSetName]["curvePtCoordsB"]
        if ptSetComm:
            nCurvePtCoordsAG = ptSetComm.allreduce(len(curvePtCoordsA), op=MPI.MAX)
            nCurvePtCoordsBG = ptSetComm.allreduce(len(curvePtCoordsB), op=MPI.MAX)
        else:
            nCurvePtCoordsAG = len(curvePtCoordsA)
            nCurvePtCoordsBG = len(curvePtCoordsB)

        # check if we need to worry about either surface
        # we will use these flags to figure out if we need to do warping.
        # we need to do the comm for the updated curves regardless
        flagA = False
        flagB = False
        if len(self.curvesOnA) > 0:
            # if len(self.surfIdxA[ptSetName]) > 0 and len(self.curvesOnA) > 0:
            flagA = True
        if len(self.curvesOnB) > 0:
            # if len(self.surfIdxB[ptSetName]) > 0 and len(self.curvesOnB) > 0:
            flagB = True

        if ptSetComm:
            rank = ptSetComm.rank
        else:
            rank = 0

        # call the bwd warping routine
        # deltaA_b is the seed for the points projected to curves
        if flagA:
            deltaA_b_local = self._warpSurfPts_b(
                dIdpt, self.points[ptSetName][0], self.surfIdxA[ptSetName], curvePtCoordsA
            )
        else:
            deltaA_b_local = numpy.zeros((N, nCurvePtCoordsAG, 3))

        # do the same for comp B
        if flagB:
            deltaB_b_local = self._warpSurfPts_b(
                dIdpt, self.points[ptSetName][0], self.surfIdxB[ptSetName], curvePtCoordsB
            )
        else:
            deltaB_b_local = numpy.zeros((N, nCurvePtCoordsBG, 3))

        # reduce seeds for both
        if ptSetComm:
            deltaA_b = ptSetComm.allreduce(deltaA_b_local, op=MPI.SUM)
            deltaB_b = ptSetComm.allreduce(deltaB_b_local, op=MPI.SUM)
        # no comm, local is global
        else:
            deltaA_b = deltaA_b_local
            deltaB_b = deltaB_b_local

        # remove the seeds for the intersection, the disps are zero (constant) so no need to diff. those
        if flagA:
            # this just returns the seeds w/o the intersection seeds
            deltaA_b = deltaA_b[:, self.nCurvePts[ptSetName]["intersection"] :]
        if flagB:
            # this just returns the seeds w/o the intersection seeds
            deltaB_b = deltaB_b[:, self.nCurvePts[ptSetName]["intersection"] :]

        # loop over the curves
        for curveName in self.featureCurveNames:

            # sizes and displacements for this curve
            sizes = self.curveProjData[ptSetName][curveName]["sizes"]
            disp = self.curveProjData[ptSetName][curveName]["disp"]

            # we get the seeds from compA seeds
            if curveName in self.curvesOnA:

                if flagA:
                    # contribution on this proc
                    deltaBar = deltaA_b[:, disp[rank] : disp[rank] + sizes[rank]].copy()

                # this proc does not have any pts projected, so set the seed to zero
                else:
                    deltaBar = numpy.zeros((N, 0, 3))

                # remove the seeds for this curve from deltaA_b seeds
                deltaA_b = deltaA_b[:, disp[-1] + sizes[-1] :]

            # seeds from compB
            elif curveName in self.curvesOnB:

                # print("[%d] checking curve: %s, flagB:"%(self.comm.rank, curveName), flagB)

                if flagB:
                    # contribution on this proc
                    deltaBar = deltaB_b[:, disp[rank] : disp[rank] + sizes[rank]].copy()

                # this proc does not have any pts projected, so set the seed to zero
                else:
                    deltaBar = numpy.zeros((N, 0, 3))

                # remove the seeds for this curve from deltaA_b seeds
                deltaB_b = deltaB_b[:, disp[-1] + sizes[-1] :]

            else:
                print("This should not happen")
            # now we have the local seeds of projection points for all functions in xyzProjb

            # get the indices of points we need to project
            idx = self.curveProjIdx[ptSetName][curveName]
            nPoints = len(idx)

            # get some data from the fwd run
            xyz = self.curveProjData[ptSetName][curveName]["xyz"]
            coor = self.curveProjData[ptSetName][curveName]["coor"]
            barsConn = self.curveProjData[ptSetName][curveName]["barsConn"]
            xyzProj = self.curveProjData[ptSetName][curveName]["xyzProj"]
            tanProj = self.curveProjData[ptSetName][curveName]["tanProj"]
            elemIDs = self.curveProjData[ptSetName][curveName]["elemIDs"]

            # we dont use tangents so seeds are zero
            tanProjb = numpy.zeros_like(tanProj)

            if nPoints > 0:
                # also get the curveMask
                curveMask = self.curveProjData[ptSetName][curveName]["curveMask"]

                # run the bwd projection for everyfunction
                for k in range(N):

                    # contribution from delta
                    xyzProjb = deltaBar[k].copy()

                    # add the contribution from dIdpt for the idx points themselves
                    xyzProjb += dIdpt[k, idx]

                    # Call fortran code (This will accumulate seeds in xyzb and self.coorb)
                    # print("[%d] right before crashing, xyzProjb.shape"%rank, xyzProjb.shape)
                    # print("[%d] right before crashing, nPoints: %d len(idx): %d, idx:"%(rank, nPoints, len(idx)), idx)
                    xyzb_new, coorb_new = curveSearchAPI.curvesearchapi.mindistancecurve_b(
                        xyz.T,
                        coor.T,
                        barsConn.T + 1,
                        xyzProj.T,
                        xyzProjb.T,
                        tanProj.T,
                        tanProjb.T,
                        elemIDs + 1,
                        curveMask,
                    )

                    # Accumulate derivatives with the correct k

                    # replace the seed in dIdpt. we subtract deltaBar here
                    # because delta was equal to the projected points minus the original points.
                    # So the delta seeds contribute with a negative sign
                    # to the seeds of the points before projection
                    dIdpt[k, idx, :] = xyzb_new.T - deltaBar[k]

                    # add the seed to the seam seed
                    self.seamBarProj[ptSetName][k, :, :] += coorb_new.T

        return compSens

    def _commCurveProj(self, pts, indices, comm):
        # this function will get the points, indices, and comm.
        # this function is called once for each feature curve
        # the indices are the indices of points that was mapped to this curve
        # we compute how many points we have mapped to this curve globally
        # furthermore, we compute the displacements
        # finally, we communicate the initial coordinates of these points,
        # that will later be used in the point-based warping

        # only do this fancy stuff if this is a "parallel" pointset
        if comm:
            nproc = comm.size

            # communicate the counts
            sizes = numpy.array(comm.allgather(len(indices)), dtype="intc")

            # total number of points
            nptsg = numpy.sum(sizes)

            # get the displacements
            disp = numpy.array([numpy.sum(sizes[:i]) for i in range(nproc)], dtype="intc")

            # sendbuf
            ptsLocal = pts[indices].flatten()
            sendbuf = [ptsLocal, len(indices) * 3]

            # recvbuf
            ptsGlobal = numpy.zeros(3 * nptsg)
            recvbuf = [ptsGlobal, sizes * 3, disp * 3, MPI.DOUBLE]

            # do an allgatherv
            comm.Allgatherv(sendbuf, recvbuf)

            # reshape into a nptsg,3 array
            curvePtCoords = ptsGlobal.reshape((nptsg, 3))

        # this is a "serial" pointset, so the results are just local
        else:
            nptsg = len(indices)
            sizes = [nptsg]
            curvePtCoords = pts[indices]

        return nptsg, sizes, curvePtCoords

    def _warpSurfPts(self, pts0, ptsNew, indices, curvePtCoords, delta):
        # this function warps points using the displacements from curve projections
        # pts0: the original surface point coordinates
        # ptsNew: updated surface pt coordinates. we will add the warped delta to these inplace
        # indices: indices of the points that we will use for this operation
        # curvePtCoords: original coordinates of points on curves
        # delta: displacements of the points on curves after projecting them

        for j in indices:

            # point coordinates with the baseline design
            # this is the point we will warp
            ptCoords = pts0[j]

            # the vectorized point-based warping we had from older versions.
            rr = ptCoords - curvePtCoords
            LdefoDist = 1.0 / numpy.sqrt(rr[:, 0] ** 2 + rr[:, 1] ** 2 + rr[:, 2] ** 2 + 1e-16)
            LdefoDist3 = LdefoDist ** 3
            Wi = LdefoDist3
            den = numpy.sum(Wi)
            interp = numpy.zeros(3)
            for iDim in range(3):
                interp[iDim] = numpy.sum(Wi * delta[:, iDim]) / den

            # finally, update the coord in place
            ptsNew[j] = ptsNew[j] + interp

    def _warpSurfPts_b(self, dIdPt, pts0, indices, curvePtCoords):

        # seeds for delta
        deltaBar = numpy.zeros((dIdPt.shape[0], curvePtCoords.shape[0], 3))

        for k in range(dIdPt.shape[0]):

            for j in indices:

                # point coordinates with the baseline design
                # this is the point we will warp
                ptCoords = pts0[j]

                # local seed for 3 coords
                localVal = dIdPt[k, j]

                # the vectorized point-based warping we had from older versions.
                rr = ptCoords - curvePtCoords
                LdefoDist = 1.0 / numpy.sqrt(rr[:, 0] ** 2 + rr[:, 1] ** 2 + rr[:, 2] ** 2 + 1e-16)
                LdefoDist3 = LdefoDist ** 3
                Wi = LdefoDist3
                den = numpy.sum(Wi)
                # interp = numpy.zeros(3)

                for iDim in range(3):
                    deltaBar[k, :, iDim] += Wi * localVal[iDim] / den

        # return the seeds for the delta vector
        return deltaBar

    def _projectToComponent(self, pts, comp, projDict):
        # we build an ADT using this component
        # from ney's code:
        # Set bounding box for new tree
        BBox = numpy.zeros((2, 3))
        useBBox = False

        # dummy connectivity data for quad elements since we have all tris
        quadConn = numpy.zeros((0, 4))

        # Compute set of nodal normals by taking the average normal of all
        # elements surrounding the node. This allows the meshing algorithms,
        # for instance, to march in an average direction near kinks.
        nodal_normals = adtAPI.adtapi.adtcomputenodalnormals(comp.nodes.T, comp.triConn.T, quadConn.T)
        comp.nodal_normals = nodal_normals.T

        # Create new tree (the tree itself is stored in Fortran level)
        adtAPI.adtapi.adtbuildsurfaceadt(
            comp.nodes.T, comp.triConn.T, quadConn.T, BBox.T, useBBox, MPI.COMM_SELF.py2f(), comp.name
        )

        # project
        numPts = pts.shape[0]
        dist2 = numpy.ones(numPts) * 1e10
        xyzProj = numpy.zeros((numPts, 3))
        normProjNotNorm = numpy.zeros((numPts, 3))

        # print("[%d] projecting to comp %s, pts.shape = "%(self.comm.rank, comp.name), pts.shape)

        # Call projection function
        procID, elementType, elementID, uvw = adtAPI.adtapi.adtmindistancesearch(
            pts.T, comp.name, dist2, xyzProj.T, comp.nodal_normals.T, normProjNotNorm.T
        )

        # Adjust indices and ordering
        elementID = elementID - 1
        uvw = uvw.T

        # normalize the normals
        normProj = tsurf_tools.normalize(normProjNotNorm)

        # deallocate ADT
        adtAPI.adtapi.adtdeallocateadts(comp.name)

        # save the data
        projDict["procID"] = procID.copy()
        projDict["elementType"] = elementType.copy()
        projDict["elementID"] = elementID.copy()
        projDict["uvw"] = uvw.copy()
        projDict["dist2"] = dist2.copy()
        projDict["normProjNotNorm"] = normProjNotNorm.copy()
        projDict["normProj"] = normProj.copy()

        # also save the original and projected points
        projDict["xyz"] = pts.copy()
        projDict["xyzProj"] = xyzProj.copy()

        # return projected points
        return xyzProj

    def _projectToComponent_b(self, dIdpt, comp, projDict):

        # we build an ADT using this component
        # from ney's code:
        # Set bounding box for new tree
        BBox = numpy.zeros((2, 3))
        useBBox = False

        # dummy connectivity data for quad elements since we have all tris
        quadConn = numpy.zeros((0, 4))

        # Compute set of nodal normals by taking the average normal of all
        # elements surrounding the node. This allows the meshing algorithms,
        # for instance, to march in an average direction near kinks.
        nodal_normals = adtAPI.adtapi.adtcomputenodalnormals(comp.nodes.T, comp.triConn.T, quadConn.T)
        comp.nodal_normals = nodal_normals.T

        # Create new tree (the tree itself is stored in Fortran level)
        adtAPI.adtapi.adtbuildsurfaceadt(
            comp.nodes.T, comp.triConn.T, quadConn.T, BBox.T, useBBox, MPI.COMM_SELF.py2f(), comp.name
        )

        # also extract the projection data we have from the fwd pass
        procID = projDict["procID"]
        elementType = projDict["elementType"]
        elementID = projDict["elementID"]
        uvw = projDict["uvw"]
        dist2 = projDict["dist2"]
        normProjNotNorm = projDict["normProjNotNorm"]
        # normProj = projDict["normProj"]

        # get the original and projected points too
        xyz = projDict["xyz"]
        xyzProj = projDict["xyzProj"]

        # also, we dont care about the normals, so the AD seeds for them should (?) be zero
        normProjb = numpy.zeros_like(normProjNotNorm)

        # also create the dIdtp for the triangulated surface nodes
        dIdptComp = numpy.zeros((dIdpt.shape[0], comp.nodes.shape[0], 3))

        # now propagate the ad seeds back for each function
        for i in range(dIdpt.shape[0]):
            # the derivative seeds for the projected points
            xyzProjb = dIdpt[i].copy()

            # Compute derivatives of the normalization process
            normProjNotNormb = tsurf_tools.normalize_b(normProjNotNorm, normProjb)

            # Call projection function
            # ATTENTION: The variable "xyz" here in Python corresponds to the variable "coor" in the Fortran code.
            # On the other hand, the variable "coor" here in Python corresponds to the variable "adtCoor" in Fortran.
            # I could not change this because the original ADT code already used "coor" to denote nodes that should be
            # projected.
            # print('\n\nbefore fortran')
            xyzb, coorb, nodal_normalsb = adtAPI.adtapi.adtmindistancesearch_b(
                xyz.T,
                comp.name,
                procID,
                elementType,
                elementID + 1,
                uvw.T,
                dist2,
                xyzProj.T,
                xyzProjb.T,
                comp.nodal_normals.T,
                normProjNotNorm.T,
                normProjNotNormb.T,
            )

            # Transpose results to make them consistent
            xyzb = xyzb.T
            coorb = coorb.T

            # Compute derivative seed contributions of the normal vectors
            # we dont need this...
            # deltaCoorb = adtAPI.adtapi.adtcomputenodalnormals_b(comp.nodes.T,
            # comp.triConn.T, quadConn.T,
            # comp.nodal_normals.T, nodal_normalsb)

            # Transpose Fortran results to make them consistent
            # deltaCoorb = deltaCoorb.T

            # put the reverse ad seed back into dIdpt
            dIdpt[i] = xyzb
            # also save the triangulated surface node seeds
            dIdptComp[i] = coorb

        # now we are done with the adt
        adtAPI.adtapi.adtdeallocateadts(comp.name)

        # call the total sensitivity of the component's dvgeo
        compSens = comp.DVGeo.totalSensitivity(dIdptComp, "triMesh")

        # the entries in dIdpt is replaced with AD seeds of initial points that were projected
        # we also return the total sensitivity contributions from components' triMeshes
        return dIdpt, compSens

    def _getUpdatedCoords(self):
        # this code returns the updated coordinates

        # first comp a
        self.compA.updateTriMesh()

        # then comp b
        self.compB.updateTriMesh()

        return

    def _getIntersectionSeam(self, comm, firstCall=False):
        # we can parallelize here. each proc gets one intersection, but needs re-structuring of some of the code.

        # this function computes the intersection curve, cleans up the data and splits the curve based on features or curves specified by the user.

        # create the dictionary to save all intermediate variables for reverse differentiation
        self.seamDict = {}

        # Call Ney's code with the quad information.
        dummyConn = numpy.zeros((0, 4))
        # compute the intersection curve, in the first step we just get the array sizes to hide allocatable arrays from python
        arraySizes = intersectionAPI.intersectionapi.computeintersection(
            self.compA.nodes.T,
            self.compA.triConn.T,
            dummyConn.T,
            self.compB.nodes.T,
            self.compB.triConn.T,
            dummyConn.T,
            self.distTol,
            comm.py2f(),
        )

        # Retrieve results from Fortran if we have an intersection
        if numpy.max(arraySizes[1:]) > 0:

            # Second Fortran call to retrieve data from the CGNS file.
            intersectionArrays = intersectionAPI.intersectionapi.retrievedata(*arraySizes)

            # We need to do actual copies, otherwise data will be overwritten if we compute another intersection.
            # We subtract one to make indices consistent with the Python 0-based indices.
            # We also need to transpose it since Python and Fortran use different orderings to store matrices in memory.
            # TODO Bug here?
            intNodes = numpy.array(
                intersectionArrays[0]
            ).T  # last entry is 0,0,0 for some reason, CHECK THIS! Checked, still zero for a proper intersection.
            # intNodes = numpy.array(intersectionArrays[0]).T[0:-1] # last entry is 0,0,0 for some reason, CHECK THIS! Checked, still zero for a proper intersection.
            barsConn = numpy.array(intersectionArrays[1]).T - 1
            parentTria = numpy.array(intersectionArrays[2]).T - 1

            # save these intermediate variables
            self.seamDict["barsConn"] = barsConn
            self.seamDict["parentTria"] = parentTria

        else:
            raise Error(
                "DVGeometryMulti Error: The components %s and %s do not intersect." % (self.compA.name, self.compB.name)
            )

        # Release memory used by fortran
        intersectionAPI.intersectionapi.releasememory()

        # sort the output
        newConn, newMap = tsurf_tools.FEsort(barsConn.tolist())

        # newConn might have multiple intersection curves
        if len(newConn) == 1:
            # we have a single intersection curve, just take this.
            seamConn = newConn[0].copy()
        # we have multiple intersection curves...
        else:
            if self.intDir is None:
                # we have multiple intersection curves but the user did not specify which direction to pick
                for i in range(len(newConn)):
                    curvename = "%s_%s_%d" % (self.compA.name, self.compB.name, i)
                    tecplot_interface.writeTecplotFEdata(intNodes, newConn[i], curvename, curvename)
                raise Error(
                    "more than one intersection curve between comps %s and %s\nThe curves are written as tecplot files in the current directory\n\nTry rerunning after specifying intDir option for the intersection."
                    % (self.compA.name, self.compB.name)
                )

            # the user did specify which direction to pick
            else:
                int_centers = numpy.zeros(len(newConn))
                # we will figure out the locations of these points and pick the one closer to the user picked direction
                for i in range(len(newConn)):
                    # get all the points
                    int_pts = intNodes[newConn[i]][:, 0]

                    # average the values
                    # the API uses a 1 based indexing, but here, we convert to a zero based indexing
                    int_centers[i] = numpy.average(int_pts[abs(self.intDir) - 1])

                # multiply the values with the sign of intDir
                int_centers *= numpy.sign(self.intDir)

                # get the argmax
                # int_index = numpy.argmax(int_centers)

                # this is the intersection seam!
                seamConn = newConn[i].copy()

        # Get the number of elements
        nElem = seamConn.shape[0]

        # now that we have a continuous, ordered seam connectivity in seamConn, we can try to detect features

        # we need to track the nodes that are closest to the supplied feature curves
        breakList = []
        curveBeg = {}
        curveBegCoor = {}

        # loop over the feature curves
        for curveName in self.featureCurveNames:

            # we need to initialize the dictionary here
            # to get the intermediate output from mindistancecurve call
            self.seamDict[curveName] = {}

            # if this curve is on compB, we use it to track intersection features
            if curveName in self.compB.barsConn and curveName not in self.remeshAll:

                # get the curve connectivity
                curveConn = self.compB.barsConn[curveName]

                # use Ney's fortran code to project the point on curve
                # first, we need to get a list of nodes that define the intersection
                intNodesOrd = intNodes[seamConn[:, 0]]

                # Get number of points
                nPoints = len(intNodesOrd)

                # Initialize references
                dist2 = numpy.ones(nPoints) * 1e10
                xyzProj = numpy.zeros((nPoints, 3))
                tanProj = numpy.zeros((nPoints, 3))
                elemIDs = numpy.zeros((nPoints), dtype="int32")

                # then find the closest point to the curve

                # only call the fortran code if we have at least one point
                # this is redundant but it is how its done in ney's code
                if nPoints > 0:
                    # This will modify xyzProj, tanProj, dist2, and elemIDs if we find better projections than dist2.
                    # Remember that we should adjust some indices before calling the Fortran code
                    # Remember to use [:] to don't lose the pointer (elemIDs is an input/output variable)
                    elemIDs[:] = (
                        elemIDs + 1
                    )  # (we need to do this separetely because Fortran will actively change elemIDs contents.
                    curveMask = curveSearchAPI.curvesearchapi.mindistancecurve(
                        intNodesOrd.T, self.compB.nodes.T, curveConn.T + 1, xyzProj.T, tanProj.T, dist2, elemIDs
                    )

                    # Adjust indices back to Python standards
                    elemIDs[:] = elemIDs - 1

                    self.seamDict[curveName]["curveMask"] = curveMask.copy()
                    self.seamDict[curveName]["elemIDs"] = elemIDs.copy()
                    self.seamDict[curveName]["intNodesOrd"] = intNodesOrd.copy()
                    self.seamDict[curveName]["xyzProj"] = xyzProj.copy()
                    self.seamDict[curveName]["tanProj"] = tanProj.copy()
                    self.seamDict[curveName]["dist2"] = dist2.copy()
                    self.seamDict[curveName]["projPtIndx"] = seamConn[:, 0][numpy.argmin(dist2)].copy()

                # now, find the index of the smallest distance
                breakList.append(numpy.argmin(dist2))

                # also get which element is the closest to the feature point
                curveBeg[curveName] = elemIDs[numpy.argmin(dist2)]

                # get which point on this element we projected to.
                curveBegCoor[curveName] = xyzProj[numpy.argmin(dist2)]

            else:
                # if it is not on compB, we still need to set up some variables so that we remesh the whole curve
                # set the beginning to the first element
                curveBeg[curveName] = 0

        # number of features we detected. This will be equal to the number of feature curves on compB
        nFeature = len(breakList)

        # if this is the first call,
        if firstCall:
            # we also save the initial curve with nodes and connectivities for distance calculations
            self.conn0 = seamConn.copy()
            self.nodes0 = intNodes.copy()
            self.nFeature = nFeature
        else:
            if nFeature != self.nFeature:
                raise Error("Number of features on the intersection curve has changed.")

        # first get an ordered list of the feature points
        # this is just our breakList "list"
        # featurePoints = intNodes[seamConn[breakList, 0]]

        # flip
        # we want breakList to be in increasing order...
        ii = 0
        for i in range(nFeature):
            # print(i,breakList[i], breakList[numpy.mod(i+1, nFeature)])
            # we loop over the breaklist elements and check if the element index is going up or down
            if breakList[i] < breakList[numpy.mod(i + 1, nFeature)]:
                ii += 1

        # print('increase index',ii)
        # now check if we need to flip the curve
        if ii == 1:  # we need at least 2 features where the element number increases...
            # we need to reverse the order of our feature curves
            # and we will flip the elements too so keep track of this change
            breakList = numpy.mod(seamConn.shape[0] - numpy.array(breakList), seamConn.shape[0])
            # TODO we had a bug in this line
            # breakList = numpy.mod(seamConn.shape[0] - numpy.flip(breakList, axis=0), seamConn.shape[0])
            # print(breakList)

            # and we need to invert the curves themselves
            seamConn = numpy.flip(seamConn, axis=0)
            seamConn = numpy.flip(seamConn, axis=1)

        # roll so that the first breakList entry is the first node
        seamConn = numpy.roll(seamConn, -breakList[0], axis=0)

        # also adjust the feature indices
        breakList = numpy.mod(breakList - breakList[0], nElem)

        # do we need this?
        # TODO figure this out?
        # print(breakList)
        # breakList.sort()

        # get the number of elements between each feature
        curveSizes = []
        for i in range(nFeature - 1):
            curveSizes.append(numpy.mod(breakList[i + 1] - breakList[i], nElem))
        # check the last curve outside the loop
        curveSizes.append(numpy.mod(breakList[0] - breakList[-1], nElem))
        # print('breaklist', breakList, "Curvesizes", curveSizes)

        # print("first node on the int:", seamConn[breakList[0], 0])
        # print("last node on the int:", seamConn[breakList[-1]+curveSizes[-1]-1, 1])

        # copy the curveSizes for the first call
        if firstCall:
            self.nElems = curveSizes[:]

        # now loop over the curves between the feature nodes. We will remesh them separately to retain resolution between curve features, and just append the results since the features are already ordered
        curInd = 0  # curveSizes[0] #+curveSizes[1]
        seam = numpy.zeros((0, 3))
        finalConn = numpy.zeros((0, 2), dtype="int32")
        for i in range(nFeature):
            # just use the same number of points *2 for now
            nNewNodes = self.nElems[i] + 1
            coor = intNodes
            barsConn = seamConn[curInd : curInd + curveSizes[i]]
            # print('i feature: %d nNewNodes: %d, curind: %d, curveSizes[i]: %d'%(i, nNewNodes, curInd, curveSizes[i]), barsConn)
            curInd += curveSizes[i]
            method = "linear"
            spacing = "linear"
            initialSpacing = 0.1
            finalSpacing = 0.1

            # re-sample the curve (try linear for now), to get N number of nodes on it spaced linearly
            # Call Fortran code. Remember to adjust transposes and indices
            newCoor, newBarsConn = utilitiesAPI.utilitiesapi.remesh(
                nNewNodes, coor.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
            )
            newCoor = newCoor.T
            newBarsConn = newBarsConn.T - 1

            # print('number of new nodes',len(newCoor))

            # add these n -resampled nodes back to back in seam and return a copy of the array
            # we don't need the connectivity info for now? we just need the coords
            # first increment the new connectivity by number of coordinates already in seam
            newBarsConn += len(seam)

            # now stack the nodes
            seam = numpy.vstack((seam, newCoor))

            # and the conn
            finalConn = numpy.vstack((finalConn, newBarsConn))

        if firstCall:
            # save the beginning and end indices of these elements
            self.seamBeg["intersection"] = 0
            self.seamEnd["intersection"] = len(finalConn)

        # save stuff to the dictionary for sensitivity computations...
        self.seamDict["intNodes"] = intNodes.copy()
        self.seamDict["seamConn"] = seamConn.copy()
        self.seamDict["curveSizes"] = curveSizes.copy()
        # size of the intersection seam w/o any feature curves
        self.seamDict["seamSize"] = len(seam)
        self.seamDict["curveBegCoor"] = curveBegCoor.copy()

        # save the intersection curve for the paper
        # if self.comm.rank == 0:
        #     curvename = '%s_%s_%d'%(self.compA.name, self.compB.name, self.counter)
        #     tecplot_interface.writeTecplotFEdata(intNodes,seamConn,curvename,curvename)

        # we need to re-mesh feature curves if the user wants...
        if self.incCurves:

            # we need to set up some variables
            if firstCall:
                self.nNodeFeature = {}
                self.distFeature = {}

            remeshedCurves = numpy.zeros((0, 3))
            remeshedCurveConn = numpy.zeros((0, 2), dtype="int32")

            # loop over each curve, figure out what nodes get re-meshed, re-mesh, and append to seam...
            for curveName in self.featureCurveNames:

                # figure out which comp owns this curve...
                if curveName in self.compB.barsConn:
                    curveComp = self.compB
                    dStarComp = self.dStarB
                elif curveName in self.compA.barsConn:
                    curveComp = self.compA
                    dStarComp = self.dStarA

                # connectivity for this curve.
                curveConn = curveComp.barsConn[curveName]

                # we already have the element that is closest to the intersection
                # or if this curve does not start from the intersection,
                # this is simply the first element
                elemBeg = curveBeg[curveName]

                # now lets split this element so that we get a better initial point...
                # this has to be on compB
                if curveName in curveBegCoor:
                    # save the original coordinate of the first point
                    ptBegSave = self.compB.nodes[curveConn[elemBeg, 0]].copy()
                    # and replace this with the starting point we want
                    self.compB.nodes[curveConn[elemBeg, 0]] = curveBegCoor[curveName].copy()

                # compute the element lengths starting from elemBeg
                firstNodes = curveComp.nodes[curveConn[elemBeg:, 0]]
                secondNodes = curveComp.nodes[curveConn[elemBeg:, 1]]
                diff = secondNodes - firstNodes
                dist2 = diff[:, 0] ** 2 + diff[:, 1] ** 2 + diff[:, 2] ** 2
                elemDist = numpy.sqrt(dist2)

                # get the cumulative distance
                cumDist = numpy.cumsum(elemDist)

                # if no marchdir for this curve, we use all of it
                if curveName in self.remeshAll:
                    # we remesh all of this curve
                    elemBeg = 0
                    elemEnd = len(curveConn)
                    if firstCall:
                        self.nNodeFeature[curveName] = elemEnd + 1
                else:
                    # do the regular thing
                    if firstCall:
                        # compute the distances from curve nodes to intersection seam
                        curvePts = curveComp.nodes[curveConn[elemBeg:, 0]]
                        # print(curvePts)

                        # Get number of points
                        nPoints = len(curvePts)

                        # Initialize references if user provided none
                        dist2 = numpy.ones(nPoints) * 1e10
                        xyzProj = numpy.zeros((nPoints, 3))
                        tanProj = numpy.zeros((nPoints, 3))
                        elemIDs = numpy.zeros((nPoints), dtype="int32")

                        # then find the closest point to the curve

                        # only call the fortran code if we have at least one point
                        if nPoints > 0:
                            # Call fortran code
                            # This will modify xyzProj, tanProj, dist2, and elemIDs if we find better projections than dist2.
                            # Remember that we should adjust some indices before calling the Fortran code
                            # Remember to use [:] to don't lose the pointer (elemIDs is an input/output variable)
                            elemIDs[:] = (
                                elemIDs + 1
                            )  # (we need to do this separetely because Fortran will actively change elemIDs contents.
                            curveMask = curveSearchAPI.curvesearchapi.mindistancecurve(
                                curvePts.T, self.nodes0.T, self.conn0.T + 1, xyzProj.T, tanProj.T, dist2, elemIDs
                            )

                        dNodes = numpy.sqrt(dist2)

                        # number of elements to use, subtract one to get the correct element count
                        nElem = (numpy.abs(dNodes - dStarComp * 1.3)).argmin() - 1

                        # we want to be one after the actual distance, so correct if needed
                        if dNodes[nElem] < dStarComp * 1.3:
                            nElem += 1

                        elemEnd = elemBeg + nElem

                        # get the total curve distance from elemBeg to this element.
                        distCurve = cumDist[nElem]

                        # save this distance as the remesh distance
                        self.distFeature[curveName] = distCurve

                        # also save how many nodes we have, we want 2 times this when re-meshing
                        self.nNodeFeature[curveName] = nElem + 1

                    else:
                        # figure out how many elements we need to go in this direction
                        elemEnd = (numpy.abs(cumDist - self.distFeature[curveName])).argmin() + elemBeg

                # get the new connectivity data between the initial and final elements
                curveConnTrim = curveConn[elemBeg : elemEnd + 1]

                # remesh the new connectivity curve, using nNode*2 times nodes
                nNewNodes = self.nNodeFeature[curveName]
                coor = curveComp.nodes
                barsConn = curveConnTrim
                method = "linear"
                spacing = "linear"
                initialSpacing = 0.1
                finalSpacing = 0.1

                # print("remeshing curve %s with barsconn:"%curveName, barsConn)
                # now re-sample the curve (try linear for now), to get N number of nodes on it spaced linearly
                # Call Fortran code. Remember to adjust transposes and indices
                newCoor, newBarsConn = utilitiesAPI.utilitiesapi.remesh(
                    nNewNodes, coor.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
                )
                newCoor = newCoor.T
                newBarsConn = newBarsConn.T - 1

                # increment the connectivitiy data
                newBarsConn += len(remeshedCurves)

                # append this new curve to the featureCurve data
                remeshedCurves = numpy.vstack((remeshedCurves, newCoor))

                remeshedCurveConn = numpy.vstack((remeshedCurveConn, newBarsConn))

                # number of new nodes added in the opposite direction
                nNewNodesReverse = 0
                if elemBeg > 0 and self.remeshBwd:
                    # also re-mesh the initial part of the curve, to prevent any negative volumes there
                    curveConnTrim = curveConn[:elemBeg]

                    nNewNodesReverse = self.nNodeFeature[curveName]
                    coor = self.compB.nodes
                    barsConn = curveConnTrim
                    method = "linear"
                    spacing = "linear"
                    initialSpacing = 0.1
                    finalSpacing = 0.1

                    # now re-sample the curve (try linear for now), to get N number of nodes on it spaced linearly
                    # Call Fortran code. Remember to adjust transposes and indices
                    newCoor, newBarsConn = utilitiesAPI.utilitiesapi.remesh(
                        nNewNodesReverse, coor.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
                    )
                    newCoor = newCoor.T
                    newBarsConn = newBarsConn.T - 1

                    newBarsConn = newBarsConn + len(remeshedCurves)

                    remeshedCurves = numpy.vstack((remeshedCurves, newCoor))
                    remeshedCurveConn = numpy.vstack((remeshedCurveConn, newBarsConn))

                if curveName in curveBegCoor:
                    # finally, put the modified initial and final points back in place.
                    # print('putting it back in fwd pass')
                    self.compB.nodes[curveConn[elemBeg, 0]] = ptBegSave.copy()

                # save some info for gradient computations later on
                self.seamDict[curveName]["nNewNodes"] = nNewNodes
                self.seamDict[curveName]["nNewNodesReverse"] = nNewNodesReverse
                self.seamDict[curveName]["elemBeg"] = elemBeg
                self.seamDict[curveName]["elemEnd"] = elemEnd
                # this includes the initial coordinates of the points for each curve
                # that has a modified initial curve
                self.seamDict["curveBegCoor"] = curveBegCoor

                if firstCall:
                    # save the beginning and end indices of these elements
                    self.seamBeg[curveName] = (
                        len(finalConn) + len(remeshedCurveConn) - (nNewNodes + nNewNodesReverse) + 2
                    )
                    self.seamEnd[curveName] = len(finalConn) + len(remeshedCurveConn)

            # now save the feature curves
            # if self.comm.rank == 0:
            #     curvename = 'featureCurves_%d'%(self.counter)
            #     tecplot_interface.writeTecplotFEdata(remeshedCurves,remeshedCurveConn,curvename,curvename)

            # now we are done going over curves,
            # so we can append all the new curves to the "seam",
            # which now contains the intersection, and re-meshed feature curves

            # increment the conn from curves
            remeshedCurveConn += len(seam)
            # stack the nodes
            seam = numpy.vstack((seam, remeshedCurves))
            # stack the conn
            finalConn = numpy.vstack((finalConn, remeshedCurveConn))

        # save the connectivity
        self.seamConn = finalConn

        # write to file to check
        # tecplot_interface.writeTecplotFEdata(seam,finalConn, 'finalcurves', 'finalcurves')

        self.counter += 1

        return seam.copy()

    def _getIntersectionSeam_b(self, seamBar, comm):
        # seamBar contains all the bwd seeds for all coordinates in self.seam

        # seam bar has shape [N, nSeamPt, 3]
        # seeds for N functions
        N = seamBar.shape[0]
        # n points in total in the combined seam
        # 3 coordinates for each point

        # allocate the space for component coordinate seeds
        coorAb = numpy.zeros((N, self.compA.nodes.shape[0], self.compA.nodes.shape[1]))
        coorBb = numpy.zeros((N, self.compB.nodes.shape[0], self.compB.nodes.shape[1]))

        # first, extract the actual intersection coordinates from feature curves.
        # we might not have any feature curves but the intersection curve will be there
        seamSize = self.seamDict["seamSize"]

        # coordinates of the beginning points for feature curves on compB
        curveBegCoor = self.seamDict["curveBegCoor"]

        intBar = seamBar[:, :seamSize, :]
        curveBar = seamBar[:, seamSize:, :]

        # dictionary to save the accumulation of curve projection seeds
        curveProjb = {}

        # check if we included feature curves
        if self.incCurves:

            # offset for the derivative seeds for this curve
            iBeg = 0

            # loop over each curve
            for curveName in self.featureCurveNames:
                # get the fwd data
                curveDict = self.seamDict[curveName]
                nNewNodes = curveDict["nNewNodes"]
                elemBeg = curveDict["elemBeg"]
                elemEnd = curveDict["elemEnd"]

                # get the derivative seeds
                newCoorb = curveBar[:, iBeg : iBeg + nNewNodes, :].copy()
                # print('newCoorb',curveName, numpy.linalg.norm(newCoorb), newCoorb)
                iBeg += nNewNodes

                # figure out which comp owns this curve...
                if curveName in self.compB.barsConn:
                    curveComp = self.compB
                    coorb = coorBb
                    # dStarComp = self.dStarB

                elif curveName in self.compA.barsConn:
                    curveComp = self.compA
                    coorb = coorAb
                    # dStarComp = self.dStarA

                # connectivity for this curve.
                curveConn = curveComp.barsConn[curveName]

                # adjust the first coordinate of the curve
                if curveName in curveBegCoor:
                    # save the original coordinate of the first point
                    ptBegSave = self.compB.nodes[curveConn[elemBeg, 0]].copy()
                    # and replace this with the starting point we want
                    self.compB.nodes[curveConn[elemBeg, 0]] = curveBegCoor[curveName].copy()

                # get the coordinates of points
                coor = curveComp.nodes

                # connectivity for this curve.
                curveConn = curveComp.barsConn[curveName]

                # trim the connectivity data
                barsConn = curveComp.barsConn[curveName][elemBeg : elemEnd + 1]

                # constant inputs
                method = "linear"
                spacing = "linear"
                initialSpacing = 0.1
                finalSpacing = 0.1

                cb = numpy.zeros((N, coor.shape[0], coor.shape[1]))

                # loop over functions
                for ii in range(N):
                    # Call Fortran code. Remember to adjust transposes and indices
                    _, _, cbi = utilitiesAPI.utilitiesapi.remesh_b(
                        nNewNodes - 1,
                        coor.T,
                        newCoorb[ii].T,
                        barsConn.T + 1,
                        method,
                        spacing,
                        initialSpacing,
                        finalSpacing,
                    )
                    # derivative seeds for the coordinates.
                    cb[ii] = cbi.T.copy()

                # print('cb',curveName, numpy.linalg.norm(cb), cb)

                # check if we adjusted the initial coordinate of the curve w/ a seam coordinate
                if elemBeg > 0:
                    if self.remeshBwd:
                        # first, we need to do the re-meshing of the other direction

                        # get the fwd data
                        nNewNodes = curveDict["nNewNodesReverse"]

                        # get the derivative seeds
                        newCoorb = curveBar[:, iBeg : iBeg + nNewNodes, :].copy()
                        # print('newCoorb',curveName, numpy.linalg.norm(newCoorb), newCoorb)
                        iBeg += nNewNodes

                        # bars conn is everything up to elemBeg
                        barsConn = curveComp.barsConn[curveName][:elemBeg]

                        # loop over functions
                        for ii in range(N):
                            # Call Fortran code. Remember to adjust transposes and indices
                            _, _, cbi = utilitiesAPI.utilitiesapi.remesh_b(
                                nNewNodes - 1,
                                coor.T,
                                newCoorb[ii].T,
                                barsConn.T + 1,
                                method,
                                spacing,
                                initialSpacing,
                                finalSpacing,
                            )
                            # derivative seeds for the coordinates.
                            cb[ii] += cbi.T.copy()

                    # the first seed is for the projected point...
                    projb = cb[:, curveConn[elemBeg, 0], :].copy()
                    # projb = cb[:,0:1,:].copy()
                    # print(curveName, projb)

                    # zero out the seed of the replaced node
                    cb[:, curveConn[elemBeg, 0], :] = numpy.zeros((N, 3))

                    # put the modified initial and final points back in place.
                    self.compB.nodes[curveConn[elemBeg, 0]] = ptBegSave.copy()

                    # we need to call the curve projection routine to propagate the seed...
                    intNodesOrd = self.seamDict[curveName]["intNodesOrd"]
                    curveMask = self.seamDict[curveName]["curveMask"]
                    elemIDs = self.seamDict[curveName]["elemIDs"]
                    xyzProj = self.seamDict[curveName]["xyzProj"]
                    tanProj = self.seamDict[curveName]["tanProj"]
                    dist2 = self.seamDict[curveName]["dist2"]

                    # we need the full bars conn for this
                    barsConn = curveComp.barsConn[curveName]

                    # allocate zero seeds
                    xyzProjb = numpy.zeros_like(xyzProj)
                    tanProjb = numpy.zeros_like(tanProj)

                    curveProjb[curveName] = numpy.zeros((N, 3))

                    for ii in range(N):

                        # the only nonzero seed is indexed by argmin dist2
                        xyzProjb[numpy.argmin(dist2)] = projb[ii].copy()

                        xyzb_new, coorb_new = curveSearchAPI.curvesearchapi.mindistancecurve_b(
                            intNodesOrd.T,
                            self.compB.nodes.T,
                            barsConn.T + 1,
                            xyzProj.T,
                            xyzProjb.T,
                            tanProj.T,
                            tanProjb.T,
                            elemIDs + 1,
                            curveMask,
                        )

                        # add the coorb_new to coorBb[ii] since coorb_new has the seeds from mindistancecurve_b
                        coorBb[ii] += coorb_new.T

                        # xyzb_new is the seed for the intersection seam node
                        # instead of saving the array full of zeros, we just save the entry we know is nonzero
                        curveProjb[curveName][ii] = xyzb_new.T[numpy.argmin(dist2)]

                # all the remaining seeds in coorb live on the component tri-mesh...
                coorb += cb

        # now we only have the intersection seam...
        nFeature = len(self.nElems)
        intNodes = self.seamDict["intNodes"]
        seamConn = self.seamDict["seamConn"]
        curveSizes = self.seamDict["curveSizes"]

        # seeds for the original intersection
        intNodesb = numpy.zeros((N, intNodes.shape[0], intNodes.shape[1]))

        # loop over each feature and propagate the sensitivities
        curInd = 0  # curveSizes[0] +curveSizes[1]
        curSeed = 0
        for i in range(nFeature):
            # just use the same number of points *2 for now
            nNewElems = self.nElems[i]
            nNewNodes = nNewElems + 1
            coor = intNodes
            barsConn = seamConn[curInd : curInd + curveSizes[i]]
            # print('SENS i feature: %d nNewNodes: %d, curind: %d, curveSizes[i]: %d'%(i, nNewNodes, curInd, curveSizes[i]))

            curInd += curveSizes[i]
            method = "linear"
            spacing = "linear"
            initialSpacing = 0.1
            finalSpacing = 0.1

            for ii in range(N):
                newCoorb = intBar[ii, curSeed : curSeed + nNewNodes, :]
                # re-sample the curve (try linear for now), to get N number of nodes on it spaced linearly
                # Call Fortran code. Remember to adjust transposes and indices
                newCoor, newBarsConn, cb = utilitiesAPI.utilitiesapi.remesh_b(
                    nNewElems, coor.T, newCoorb.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
                )
                intNodesb[ii] += cb.T
                # print('i=%d ii=%d, len(cb)'%(i,ii),len(cb.T))
                # print(cb.T)

            # test if newCoor is what we actually have in the intersection
            # print("sizes of the seed and pts")
            # print('i feature: %d'%i)
            # print("shape of intBar",intBar.shape)
            # print('curSeed : curseed+nNewNodes', curSeed, curSeed+nNewNodes)
            # print('checking newcoor')
            # print(newCoor.T - self.seam[curSeed:curSeed+nNewNodes])

            curSeed += nNewNodes

        # add the contributions from the curve projection if we have any
        for curveName, v in curveProjb.items():
            # get the index
            idx = self.seamDict[curveName]["projPtIndx"]
            # print('adding contributions from', curveName, v)
            for ii in range(N):
                # add the contribution
                intNodesb[ii, idx] += v[ii]

        dummyConn = numpy.zeros((0, 4))
        barsConn = self.seamDict["barsConn"]
        parentTria = self.seamDict["parentTria"]
        # do the reverse intersection computation to get the seeds of coordinates
        for ii in range(N):
            cAb, cBb = intersectionAPI.intersectionapi.computeintersection_b(
                self.compA.nodes.T,
                self.compA.triConn.T,
                dummyConn.T,
                self.compB.nodes.T,
                self.compB.triConn.T,
                dummyConn.T,
                intNodes.T,
                intNodesb[ii].T,
                barsConn.T + 1,
                parentTria.T + 1,
                self.distTol,
            )

            coorAb[ii] += cAb.T
            coorBb[ii] += cBb.T

        # get the total sensitivities from both components
        compSens_local = {}
        compSensA = self.compA.DVGeo.totalSensitivity(coorAb, "triMesh")
        for k, v in compSensA.items():
            kNew = "%s:%s" % (self.compA.name, k)
            compSens_local[kNew] = v

        compSensB = self.compB.DVGeo.totalSensitivity(coorBb, "triMesh")
        for k, v in compSensB.items():
            kNew = "%s:%s" % (self.compB.name, k)
            compSens_local[kNew] = v

        # finally sum the results across procs if we are provided with a comm
        if comm:
            compSens = {}
            # because the results are in a dictionary, we need to loop over the items and sum
            for k, v in compSens_local.items():
                compSens[k] = comm.allreduce(v, op=MPI.SUM)
        else:
            # we can just pass the dictionary
            compSens = compSens_local

        return compSens
