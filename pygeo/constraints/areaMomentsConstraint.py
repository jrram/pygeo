# ======================================================================
#         Imports
# ======================================================================
import numpy as np
from baseConstraint import GeometricConstraint # TODO: had to remove dot for this to work; change it back


class AreaMomentsConstraint(GeometricConstraint):
    """
    This class is used to represet individual second moment of area constraint.
    """
    # The parameter list is explained in the addVolumeConstaint() of
    # the DVConstraints class
    # """

    def __init__(self, name, nSpan, nChord, coords, lower, upper, scaled, scale, DVGeo, addToPyOpt):

        self.name = name
        self.nCon = 1
        self.nSpan = nSpan
        self.nChord = nChord
        self.coords = coords
        self.lower = lower
        self.upper = upper
        self.scaled = scaled
        self.scale = scale
        self.DVGeo = DVGeo
        self.addToPyOpt = addToPyOpt
        self.flipVolume = False

        GeometricConstraint.__init__(
            self, self.name, self.nCon, self.lower, self.upper, self.scale, self.DVGeo, self.addToPyOpt
        )

        # First thing we can do is embed the coordinates into DVGeo
        # with the name provided:
        self.DVGeo.addPointSet(self.coords, self.name)

        # Now get the reference second moments of area
        self.IJ0 = self.evalAreaMoments()

# IJ should have three components: Ixx, Iyy, and J

    def evalFunctions(self, funcs, config):
        """
        Evaluate the function this object has and place in the funcs dictionary

        Parameters
        ----------
        funcs : dict
            Dictionary to place function values
        """
        # Pull out the most recent set of coordinates:
        self.coords = self.DVGeo.update(self.name, config=config)
        Ixx, Iyy, Jz = self.evalAreaMoments()
        if self.scaled:
            IJ /= self.IJ0 # TODO: check division
        funcs[f"{self.name}_Ixx"] = Ixx
        funcs[f"{self.name}_Iyy"] = Iyy
        funcs[f"{self.name}_Jz"] = Jz

    def evalFunctionsSens(self, funcsSens, config):
        """
        Evaluate the sensitivity of the functions this object has and
        place in the funcsSens dictionary

        Parameters
        ----------
        funcsSens : dict
            Dictionary to place function values
        """
        pass

    def evalAreaMoments(self):
        """
        Evaluate the second moments of area for each section
        """
        #IJ = {} # TODO: could make these three into a dictionary of arrays with Ixx, etc key?
        Ixx = np.zeros(self.nSpan - 1) # TODO check that size is correct
        Iyy = np.zeros(self.nSpan - 1)
        Jz = np.zeros(self.nSpan - 1) 
        x = self.coords.reshape((self.nSpan, self.nChord, 2, 3))
        for i in range(self.nSpan - 1):

            # Compute centroid of airfoil
            xa = 0
            ya = 0
            for j in range(self.nChord - 1):
                # Add centroid of each strip
                xa += 0.5 * (x[i, j+1, 0] + x[i, j, 0])
                ya += 0.5 * (x[i, j, 1] + x[i, j, 0])
            # Divide totals by number of strips to get average
            xa /= self.nChord - 1 # TODO: check number
            ya /= self.nChord - 1 # TODO: check number

            # Compute second moments of area of airfoil  
            for j in range(self.nChord - 1):

                # Calculate width, height and centroid of current strip
                dx = x[i, j+1, 0] - x[i, j, 0]
                dy = x[i, j, 1] - x[i, j, 0]
                xc = 0.5 * (x[i, j+1, 0] + x[i, j, 0]) # TODO: store these from previous xc, yc computation?
                yc = 0.5 * (x[i, j, 1] + x[i, j, 0])

                # Debuggin; Checks prints
                # Are dx and dy equal on both sides as assumed?
                print("dx:", dx, "=?", x[i, j+1, 1] - x[i, j, 1])
                print("dy:", yx, "=?", x[i, j+1, 1] - x[i, j+1, 0])

                # Debugging: Write data for each strip

                # Calculate second moment of area of each strip
                Ixx[i] += (dx * dy**3) / 12 + (yc - ya)**2 * dx * dy
                Iyy[i] += (dx**3 * dy) / 12 + (xc - xa)**2 * dx * dy
            # Calculate polar second moment of area
            Jz[i] = Ixx + Iyy

        return Ixx, Iyy, Jz

    # def evalVolumeSens(self):
    #     """
    #     Evaluate the derivative of the volume with respect to the
    #     coordinates
    #     """

    # def evalFunctionsSens(self, funcsSens, config):
    #     """
    #     Evaluate the sensitivity of the functions this object has and
    #     place in the funcsSens dictionary

    #     Assumes that evalFunctions method was called just prior
    #     and results stashed on proc 0

    #     Parameters
    #     ----------
    #     funcsSens : dict
    #         Dictionary to place function values
    #     config : str
    #         The DVGeo configuration to apply this constraint to. Must be either None
    #         which will apply to *ALL* the local DV groups or a single string specifying
    #         a particular configuration.
    #     """


    
    def writeTecplot(self, handle):
        """No need to write the composite volume since each of the
        individual ones are already written"""
        pass

class dummyDVgeo:

    def addPointSet(self, *args):
        pass

if __name__ == "__main__":
    print("Hello!")
    nSpan = 1
    nChord = 2
    x = np.zeros([nSpan, nChord, 2, 3]) # nSpan, nChord, top/botom, xyz indices (x=0, y=1)
    x[0,0,0,:] = np.array([0,0,0]) # bottom left corner
    x[0,1,0,:] = np.array([1,0,0]) # bottom right corner
    x[0,0,1,:] = np.array([0,1,0]) # top left corner
    x[0,1,1,:] = np.array([1,1,0]) # top right corner
    coords = x.reshape((nSpan * nChord * 2, 3))
    lower = -1.e20
    upper = 1.e20
    scaled = False
    scale = 1.0
    DVGeo = dummyDVgeo()
    addToPyOpt = True
    foil = AreaMomentsConstraint("foil", 1, 2, coords, lower, upper, scaled, scale, DVGeo, addToPyOpt)