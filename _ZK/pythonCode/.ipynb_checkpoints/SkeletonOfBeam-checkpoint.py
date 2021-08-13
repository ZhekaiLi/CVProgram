import numpy as np
import sympy as sy
import scipy as sp

class SkeletonOfBeam:
    """
    Get skeleton of beam

    Arributes:
        mesh: Input trimesh.mesh
        centroid: An array of the coordinates of mesh's centroid
        nVec: Initial skeleton vector
        
        Intersections: A list of trimesh.Path3D
        SkeletonPoints: A list of np.array
        
        XYZCoordinate: A new coordinate which use nVec as the x-axis
        XYProjections: projections of SkeletonPoints on the x-y plane of XYZCoordinate
        XZProjections: projections of SkeletonPoints on the x-z plane of XYZCoordinate
        
    Used packages:
        import numpy as np
        import sympy as sy
        import scipy as sp
    """
    mesh = None
    centroid = None
    nVec = None
    
    IntersectionScale = None # scale along the nVec to get Intersections
    Intersections = None
    Intersections2D = None
    SkeletonPoints = None # centroids of the intersections
    
    XYZCoordinate = None # 3*3 array, first line = x-axis ...
    XYProjections = None
    XZProjections = None
    
    u_xyPlane = None
    u_xzPlane = None
    alpha_xyPlane = None
    alpha_xzPlane = None
    # Sympy derivation of the projection of skeleton curve on x-y plane
    # Calculate the derivation value by using dudx_xyPlane.evalf(subs={'xi': xi_value})
    # where xi_value in [0, 1]
    dudx_xyPlane = None
    dudx_xzPlane = None # similar to dudx_xyPlane
    
    
    
    def __init__(self, mesh, rough_normalVector):
        self.mesh = mesh
        self.centroid = mesh.centroid.copy()
        # 质心截面的一个大致的法向量（目前单指 [1, 0, 0]）
        self.nVec = np.asarray(rough_normalVector)
        
    
    def getScaleAlongSkeletonVec(self):
        """Get the scale along the input skeleton vector, which then serve for taking intersections
        """
        currentScale = 0.2
        
        # Get a rough scale
        ifSliceExist = True     
        while ifSliceExist:
            originPoint = self.centroid + currentScale * self.nVec
            slice = self.mesh.section(plane_origin=originPoint,  plane_normal=self.nVec) # take the slice
            if slice is None:
                ifSliceExist = False
            else:
                currentScale = currentScale * 2
             
        # Get the accurate sacles
        scales = [0, 0]
        
        # Binary Search for positive half
        leftS = currentScale / 2
        rightS = currentScale
        while rightS - leftS > 1e-2:
            midS = (leftS + rightS) / 2
            originPoint = self.centroid + midS * self.nVec
            slice = self.mesh.section(plane_origin=originPoint,  plane_normal=self.nVec)
            if slice is None: rightS = midS
            else: leftS = midS
        scales[1] = leftS
                
        # Binary Search for negative half
        leftS = -currentScale
        rightS = -currentScale / 2
        while rightS - leftS > 1e-2:
            midS = (leftS + rightS) / 2
            originPoint = self.centroid + midS * self.nVec
            slice = self.mesh.section(plane_origin=originPoint,  plane_normal=self.nVec)
            if slice is None: leftS = midS
            else: rightS = midS
        scales[0] = rightS
        
        self.IntersectionScale = np.asarray(scales)
            
        
          
    def getIntersectionsFromStep(self, step=1):
        """Get the intersections of beam along nVec
        
        :param step: interval of intersections
        """
        sections = []
        extents = self.IntersectionScale # 截取区间
        levels = np.arange(*extents, step=step)  # 每隔 1m 截一次
        for i in range(len(levels)):
            origin_temp = self.centroid.copy()
            origin_temp[0] = origin_temp[0] + levels[i]
            try:
                slice = self.mesh.section(plane_origin=origin_temp,  plane_normal=self.nVec)
                # 选取每个截面图中面积最大的子图，实现初步去噪
                if slice is not None:
                    slice_2D, to_3D = slice.to_planar()
                    slices_splited = slice_2D.split()
                    sliceIndex = np.argmax([s.area for s in slices_splited])
                    slice_2D = slices_splited[sliceIndex]
                    if slice_2D.area > 1e-1:
                        sections.append(slice_2D.to_3D(to_3D))
            except:
                pass
        
        self.Intersections = sections
        
    def getIntersectionsFromSliceNum(self, sliceNum=5):
        """Get the intersections of beam along nVec
        
        :param sliceNum: number of intersections
        """
        sections = []
        extents = self.IntersectionScale # 截取区间
        levels = np.linspace(*extents, sliceNum)  # 每隔 1m 截一次
        for i in range(len(levels)):
            origin_temp = self.centroid.copy()
            origin_temp[0] = origin_temp[0] + levels[i]
            try:
                slice = self.mesh.section(plane_origin=origin_temp,  plane_normal=self.nVec)
                # 选取每个截面图中面积最大的子图，实现初步去噪
                if slice is not None:
                    slice_2D, to_3D = slice.to_planar()
                    slices_splited = slice_2D.split()
                    sliceIndex = np.argmax([s.area for s in slices_splited])
                    slice_2D = slices_splited[sliceIndex]
                    if slice_2D.area > 1e-1:
                        sections.append(slice_2D.to_3D(to_3D))
            except:
                pass
        
        self.Intersections = sections
    
                
    def getSkeletonPoints(self):   
        """Get the centroids of intersections
        """
        self.SkeletonPoints = []
        for s in self.Intersections:
            self.SkeletonPoints.append(s.centroid)
            
    def getNewCoordinate(self):
        """Get new XYZ coordinate with vector nVec as x-axis
        """
        R = GeometryToolBox.rotation_matrix([1, 0, 0], self.nVec)
        x = self.nVec
        y = np.dot(R, np.array([0, 1, 0]))
        z = np.dot(R, np.array([0, 0, 1]))
        self.XYZCoordinate = np.vstack((x, y, z))
        return self.XYZCoordinate
    
    def getProjections(self):
        """Get the projections of the SkeletonPoints on x-y, x-z plane
        """    
        x, y, z = self.XYZCoordinate
        origin = self.SkeletonPoints[0]
        self.XYProjections = [GeometryToolBox.projected_point(p, origin, x, y) for p in self.SkeletonPoints]
        self.XZProjections = [GeometryToolBox.projected_point(p, origin, x, z) for p in self.SkeletonPoints]
        
    
    def getSkeletonEqs(self):
        """Get the equation of the skeleton curve in x-y, x-z plane
        """
        xs = np.array(self.XYProjections)[:,0]
        ys = np.array(self.XYProjections)[:,1]
        zs = np.array(self.XZProjections)[:,1]

        L = xs[-1] - xs[0]
        xis = xs / L

        errorValue = lambda x,y,A: y - np.dot(A, x)
        a_init = np.array([1] * 4)

        # Calculate the derivation equation on x-y plane
        # Get the optimal parameters using least squre error method
        a1 = sp.optimize.leastsq(errorValue, a_init, args=(ys, self._H(xis, L)))[0]
        self.alpha_xyPlane = a1
        
        # Derivation
        xi = sy.symbols('xi')
        self.u_xyPlane = (self._H(xi, L, ifsymbol=True) * a1).sum()
        
        # Then calculate the derivation equation on x-z plane
        a2 = sp.optimize.leastsq(errorValue, a_init, args=(zs, self._H(xis, L)))[0]
        self.alpha_xzPlane = a2
        self.u_xzPlane = (self._H(xi, L, ifsymbol=True) * a2).sum()
        
    
    def getDerivativeSkeletonEqs(self):
        """Get the derivation of the skeleton curve in x-y, x-z plane
        """
        xs = np.array(self.XYProjections)[:,0]
        L = xs[-1] - xs[0]
        
        # Derivation
        xi = sy.symbols('xi')
        self.dudx_xyPlane = sy.diff(self.u_xyPlane, xi) / L
        
        # Then calculate the derivation equation on x-z plane
        self.dudx_xzPlane = sy.diff(self.u_xzPlane, xi) / L
        
    def getNewSkeletonPoints(self):
        """Get new centroids according to the curve equations
        """
        xVec, yVec, zVec = self.XYZCoordinate
        xs = np.array(self.XYProjections)[:,0]
        L = xs[-1] - xs[0]
        xis = xs / L

        self.SkeletonPoints = []
        for i in range(len(xis)):
            xi_value = xis[i]
            sX = xs[i]
            sY = self.u_xyPlane.evalf(subs={'xi': xi_value})
            sZ = self.u_xyPlane.evalf(subs={'xi': xi_value})
            self.SkeletonPoints.append(sX*xVec + sY*yVec + sZ*zVec)
        
    def getNewIntersections(self):
        """Get the intersections of beam along vector [1, 0, 0]
        
        :param step: interval of intersections
        """
        sections = []
        sections2D = []
        xs = np.array(self.XYProjections)[:,0]
        L = xs[-1] - xs[0]
        xis = xs / L
        if len(xis) != len(self.SkeletonPoints):
            raise Exception("Conflit between xis and SkeletonPoints.", self.SKeletonPoints)
        for i in range(len(xis)):
            xi = xis[i]
            normalVec = self.returnTangentVectorAtXi(xi)
            originPoint = self.SkeletonPoints[i]
            try:
                slice = self.mesh.section(plane_origin=originPoint,  plane_normal=normalVec)
                # 选取每个截面图中面积最大的子图，实现初步去噪
                if slice is not None:
                    slice_2D, to_3D = slice.to_planar()
                    slices_splited = slice_2D.split()
                    sliceIndex = np.argmax([s.area for s in slices_splited])
                    slice_2D = slices_splited[sliceIndex]
                    if slice_2D.area > 1e-1:
                        sections2D.append(slice_2D)
                        sections.append(slice_2D.to_3D(to_3D))
            except:
                pass
        
        self.Intersections2D = sections2D
        self.Intersections = sections            
        
        
    def returnTangentVectorAtXi(self, xi_value):
        """Return the tangent vector at xi=xi_value
        """
        xVec, yVec, zVec = self.XYZCoordinate

        sY = self.dudx_xyPlane.evalf(subs={'xi': xi_value}) # scalar y
        sZ = self.dudx_xzPlane.evalf(subs={'xi': xi_value}) # scalar z
        
        tanVec = xVec + sY*yVec + sZ*zVec
        
        return tanVec.astype(np.float64)
    
    def returnSkeletonPointsInXiRange(self, xi_s, xi_e, pNum=10):
        points = []
        
        xVec, yVec, zVec = self.XYZCoordinate
        xs = np.array(self.XYProjections)[:,0]
        L = xs[-1] - xs[0]
        
        xis = np.linspace(xi_s, xi_e, pNum)
        xs = xis * L
        for i in range(len(xis)):
            xi_value = xis[i]
            sX = xs[i]
            sY = self.u_xyPlane.evalf(subs={'xi': xi_value})
            sZ = self.u_xyPlane.evalf(subs={'xi': xi_value})
            points.append(sX*xVec + sY*yVec + sZ*zVec)
        
        return np.asarray(points).astype(np.float64)
        
        
    
    def showIntersections(self, ifOverlapped=True):
        """Visualize the intersections in one graph
        params:
        ifOverlapped: if True, all intersections will be plotted into one graph
            else Flase, each intersection will be plotted in order
        """
        if ifOverlapped:
            combined = np.sum(self.Intersections2D)
            combined.show()
        else:
            for c2 in self.Intersections2D:
                c2.show()
                
    
    def _H(self, xs, L, ifsymbol=False): 
        h2 = 1 - 3*xs**2 + 2*xs**3
        h3 = xs*(1 - 2*xs + xs**2)*L
        h5 = xs**2*(3 - 2*xs)
        h6 = xs**2*(xs - 1)*L
        if ifsymbol:
            return np.array([h2, h3, h5, h6])
        else:
            return np.hstack((h2.reshape(len(xs), -1), h3.reshape(len(xs), -1), 
                              h5.reshape(len(xs), -1), h6.reshape(len(xs), -1)))
        
        
#----------------------------------------------------------------------------------------------------

class GeometryToolBox:
    def rotation_matrix(vS, vE):
        """
        """
        vS = np.asarray(vS) / np.linalg.norm(vS) # 归一化
        vE = np.asarray(vE) / np.linalg.norm(vE)
        axis = np.cross(vS, vE)
        theta = np.arccos(np.dot(vS, vE))

        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d

        return np.array([[aa+bb-cc-dd, 2*(bc + ad), 2*(bd - ac)],
                         [2*(bc - ad), aa+cc-bb-dd, 2*(cd + ab)],
                         [2*(bd + ac), 2*(cd - ab), aa+dd-bb-cc]])


    def projected_point(point, plane_origin, planeVec1, planeVec2):
        """
        已知平面内一原点以及两个正交向量，求已知点在该平面内的投影坐标（二维）
        """
        pVec = np.asarray(point) - np.asarray(plane_origin)
        xVec = np.asarray(planeVec1) / np.linalg.norm(planeVec1)
        yVec = np.asarray(planeVec2) / np.linalg.norm(planeVec2)
        s1 = np.dot(pVec, xVec)
        s2 = np.dot(pVec, yVec)

        return [s1, s2]