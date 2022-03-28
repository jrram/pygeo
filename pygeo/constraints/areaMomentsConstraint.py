# ======================================================================
#         Imports
# ======================================================================
import numpy as np
from .baseConstraint import GeometricConstraint

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
        self.Ixx0, self.Iyy0, self.Jz0 = self.evalAreaMoments()

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
            Ixx = np.divide(Ixx, self.Ixx0)
            Iyy = np.divide(Iyy, self.Iyy0)
            Jz = np.divide(Jz, self.Jz0)
        funcs[f"{self.name}_Ixx"] = Ixx
        funcs[f"{self.name}_Iyy"] = Iyy
        funcs[f"{self.name}_Jz"] = Jz

    # def evalFunctionsSens(self, funcsSens, config):
    #     """
    #     Evaluate the sensitivity of the functions this object has and
    #     place in the funcsSens dictionary

    #     Parameters
    #     ----------
    #     funcsSens : dict
    #         Dictionary to place function values
    #     """
    #     nDV = self.DVGeo.getNDV()
    #     if nDV > 0:
    #         dVdPt = self.evalVolumeSens()
    #         if self.scaled:
    #             dVdPt /= self.V0

    #         # Now compute the DVGeo total sensitivity:
    #         funcsSens[self.name] = self.DVGeo.totalSensitivity(dVdPt, self.name, config=config)

    def evalAreaMoments(self):
        """
        Evaluate the second moments of area for each section
        """
        #IJ = {} # TODO: could make these three into a dictionary of arrays iwth Ixx, etc lable?
        Ixx = np.zeros(self.nSpan)
        Iyy = np.zeros(self.nSpan)
        Jz = np.zeros(self.nSpan)
        # Reshape coordinates
        coords = self.coords.reshape((self.nSpan, self.nChord, 2, 3))
        # [] TODO JM-: For RAE2822 airfoil, I had to flip the top/bottom index (3rd index below) to get positive dy. Need to test again in box and other cases. Is there a consistent definition?
        # Chordwise coordinate
        x_bot = coords[:,:,1,0]
        x_top = coords[:,:,0,0]
        # Vertical coordinate
        y_bot = coords[:,:,1,1]
        y_top = coords[:,:,0,1]

        print(f"\nself.nSpan={self.nSpan}")
        print(f"self.nChord={self.nChord}")
        print(f"x_top = {x_top}")
        print(f"x_bot = {x_bot}")
        print(f"y_top = {y_top}")
        print(f"y_bot = {y_bot}")

        # Loop over each spanwise location
        for i in range(self.nSpan):

            # Airfoil centroid
            xc = 0
            yc = 0
            # Airfoil area and chord
            area = 0
            chord = 0

            # Loop over chordwise indices for current airfoil
            for j in range(self.nChord - 1):

                # Rectangular strip with midpoint values
                ## Sides
                dx = x_bot[i,j+1] - x_bot[i,j]
                chord += dx
                assert dx > 0, "dx is not positive!"
                assert dx == x_top[i,j+1] - x_top[i,j], "dx not same on top and bottom!"
                ## y coordinate in middle top of rectangle
                y_rect_top = (y_top[i,j] + y_top[i,j+1]) / 2
                # y coordinate in middle bottom of rectangle
                y_rect_bot = (y_bot[i,j] + y_bot[i,j+1]) / 2
                # Height of rectangle
                dy_rect = y_rect_top - y_rect_bot
                assert dy_rect > 0, "dy_rect is not positive!"
                # Area of rectangular strip
                area_rect = dx * dy_rect
                area += area_rect
                # Centroid of rectangular strip
                xc_rect = x_bot[i,j] + dx/2
                yc_rect = y_rect_bot + dy_rect/2
                # Add centroid of rectangular strip
                xc += xc_rect * area_rect
                yc += yc_rect * area_rect
                # Second moments of area about global coordinates
                Ixx[i] += (dx * dy_rect**3) / 12 + area_rect * yc_rect**2
                Iyy[i] += (dx**3 * dy_rect) / 12 + area_rect * xc_rect**2

                # # Drela's version of torsional rigidity
                Jz[i] += dy_rect**3 * dx / 3
                
            # End chordwise loop

            # Compute centroid
            xc /= area
            yc /= area
            # Convert moments of area to be relative to airfoil centroid
            Ixx[i] -= area * yc**2
            Iyy[i] -= area * xc**2

            # Approximate torsional rigidity based on one rectangular strip
            # t_avg = area / chord
            # Jz[i] = chord * t_avg**3 / 3

            print(f"\nArea = {area:.6e}, \nChord = {chord:.6e}")
            print(f"Centroid = {xc:.6}, {yc:.6}")
            print(f"i={i:3d}, Ixx={Ixx[i]:.6e}, Iyy={Iyy[i]:.6e}, Jz={Jz[i]:.6e}")

        # End spanwise loop
        return Ixx, Iyy, Jz

    # def evalAreaMomentsSens(self):
    #     """
    #     Evaluate the derivative of the moments of area with respect to the
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

if __name__ == "__main__":

    print("Testing area moments...")

    # Simple rectangle
    print('\nTest 1: One rectangle')
    #coords = self.coords.reshape((self.nSpan, self.nChord, 2, 3))
    coords = np.zeros([1, 2, 2, 3])
    print(coords[0,0,0,:].shape)
    coords[0,0,0,:] = [1,2,0]
    coords[0,1,0,:] = [3,2,0]
    coords[0,0,1,:] = [1,6,0]
    coords[0,1,1,:] = [3,6,0]
    acon = AreaMomentsConstraint('acon', 1,2,coords) #,0,0,False,1,None,False)
    #acon.evalFunctions()
    Ixx, Iyy, Jz = acon.evalAreaMoments()
    print(Ixx, Iyy, Jz)

    # One rectangle plus triangles
    print('\nTest 2: One rectangle plus triangles')
    coords = np.zeros([1, 2, 2, 3])
    #print(coords[0,0,0,:].shape)
    coords[0,0,0,:] = [1,2,0]
    coords[0,1,0,:] = [3,1.5,0]
    coords[0,0,1,:] = [1,6,0]
    coords[0,1,1,:] = [3,6.5,0]
    acon = AreaMomentsConstraint('acon', 1,2,coords) #,0,0,False,1,None,False)
    #acon.evalFunctions()
    Ixx, Iyy, Jz = acon.evalAreaMoments()
    print(Ixx, Iyy, Jz)

    # Two rectangles plus triangles
    print('\nTest 3: Two rectangles plus triangles')
    coords = np.zeros([1, 3, 2, 3])
    #print(coords[0,0,0,:].shape)
    coords[0,0,0,:] = [1,2,0]
    coords[0,1,0,:] = [3,1.5,0]
    coords[0,2,0,:] = [5,2,0]
    coords[0,0,1,:] = [1,6,0]
    coords[0,1,1,:] = [3,6.5,0]
    coords[0,2,1,:] = [5,6,0]
    acon = AreaMomentsConstraint('acon', 1, 3, coords) #,0,0,False,1,None,False)
    #acon.evalFunctions()
    Ixx, Iyy, Jz = acon.evalAreaMoments()
    print(Ixx, Iyy, Jz)