import numpy as np
import open3d as o3d
import matplotlib.pyplot as mp

class find_end:
    """
    To find 4 end points of the beam

    Arributes:
         mesh: input mesh---o3d.io.read_triangle_mesh("cantilever.obj")
         elevation: input angle in degree
         zimuth: input angle  in degree

         point: vertices of the mesh 
         point2: hull vertices 
         



    Used packages:
        import numpy as np
        import open3d as o3d
        import matplotlib.pyplot as mp


    """
    
    
    def __init__(self,mesh):
        self.pcd=mesh
        self.point=np.asarray(self.pcd.vertices)

    def get_hull(self):
        self.hull, _ = self.pcd.compute_convex_hull()
        self.hull_ls = o3d.geometry.LineSet.create_from_triangle_mesh(self.hull)
        self.point2=np.asarray(self.hull.vertices)
    
    def show_hull(self):
        self.hull_ls.paint_uniform_color((1, 0, 0))
        o3d.visualization.draw_geometries([self.pcd, self.hull_ls])
    
    def show_points(self):
        self.x=self.point2[:,0]
        self.y=self.point2[:,1]
        self.z=self.point2[:,2]

        # 2.绘制图片
        mp.figure("3D Scatter", facecolor="lightgray")
        ax3d = mp.gca(projection="3d")  # 创建三维坐标

        mp.title('possible points', fontsize=20)
        ax3d.set_xlabel('x', fontsize=14)
        ax3d.set_ylabel('y', fontsize=14)
        ax3d.set_zlabel('z', fontsize=14)
        mp.tick_params(labelsize=10)



        ax3d.scatter(self.x,self.y,self.z, s=20,  cmap="jet", marker="o")

        mp.show()

    def show_leftpoints(self,elevation,zimuth):
        n=self.x.size
        self.x_left=[]
        self.y_left=[]
        self.z_left=[]
        self.x_right=[]
        self.y_right=[]
        self.z_right=[]
        for i in range(n):
            if self.x[i]<-10:
                self.x_left.append(self.x[i])
                self.y_left.append(self.y[i])
                self.z_left.append(self.z[i])
            else:
                self.x_right.append(self.x[i])
                self.y_right.append(self.y[i])
                self.z_right.append(self.z[i])



        mp.figure("3D Scatter", facecolor="lightgray")
        ax3d = mp.gca(projection="3d")  # 创建三维坐标

        mp.title('left points', fontsize=20)
        ax3d.set_xlabel('x', fontsize=14)
        ax3d.set_ylabel('y', fontsize=14)
        ax3d.set_zlabel('z', fontsize=14)
        mp.tick_params(labelsize=10)



        ax3d.view_init(elev=elevation,    # 仰角
                    azim=zimuth       # 方位角
                    )


        ax3d.scatter(self.x_left,self.y_left,self.z_left, s=20,  cmap="jet", marker="o")


        mp.show()   

    def show_2dleftpoints(self):
        mp.subplot(221)
        mp.scatter(self.x_left,self.y_left, s=20,  cmap="jet", marker="o")

        mp.title("x-y") #图名
        mp.xlabel("x")#x轴标签
        mp.ylabel("y")#y轴标签
        mp.tick_params(axis='both')#x,y轴都有刻度

        mp.subplot(222)
        mp.scatter(self.x_left,self.z_left, s=20,  cmap="jet", marker="o")

        mp.title("x-z") #图名
        mp.xlabel("x")#x轴标签
        mp.ylabel("z")#y轴标签
        mp.tick_params(axis='both')#x,y轴都有刻度

        mp.subplot(223)
        mp.scatter(self.y_left,self.z_left, s=20,  cmap="jet", marker="o")

        mp.title("y-z") #图名
        mp.xlabel("y")#x轴标签
        mp.ylabel("z")#y轴标签
        mp.tick_params(axis='both')#x,y轴都有刻度

        mp.show()
    
    def show_leftend(self,xl,yl,zl):
        mp.figure("3D Scatter", facecolor="lightgray")
        ax3d = mp.gca(projection="3d")  # 创建三维坐标

        mp.title('left points', fontsize=20)
        ax3d.set_xlabel('x', fontsize=14)
        ax3d.set_ylabel('y', fontsize=14)
        ax3d.set_zlabel('z', fontsize=14)
        mp.tick_params(labelsize=10)



        ax3d.view_init(elev=0,    # 仰角
                    azim=30       # 方位角
                    )


        ax3d.scatter(self.x_left,self.y_left,self.z_left, s=20,  cmap="jet", marker="o")

        ax3d.scatter(xl,yl,zl, s=20,  cmap='red', marker="o")

        mp.show()

    def show_rightpoints(self,elevation,zimuth):
        mp.figure("3D Scatter", facecolor="lightgray")
        ax3d = mp.gca(projection="3d")  # 创建三维坐标

        mp.title('right points', fontsize=20)
        ax3d.set_xlabel('x', fontsize=14)
        ax3d.set_ylabel('y', fontsize=14)
        ax3d.set_zlabel('z', fontsize=14)
        mp.tick_params(labelsize=10)

        ax3d.view_init(elev=elevation,    # 仰角
                    azim=zimuth       # 方位角
                    )


        ax3d.scatter(self.x_right,self.y_right,self.z_right, s=20,  cmap="jet", marker="o")


        mp.show()


    def show_2drightpoints(self):
        mp.subplot(221)
        mp.scatter(self.x_right,self.y_right, s=20,  cmap="jet", marker="o")

        mp.title("x-y") #图名
        mp.xlabel("x")#x轴标签
        mp.ylabel("y")#y轴标签
        mp.tick_params(axis='both')#x,y轴都有刻度

        mp.subplot(222)
        mp.scatter(self.x_right,self.z_right, s=20,  cmap="jet", marker="o")

        mp.title("x-z") #图名
        mp.xlabel("x")#x轴标签
        mp.ylabel("z")#y轴标签
        mp.tick_params(axis='both')#x,y轴都有刻度

        mp.subplot(223)
        mp.scatter(self.y_right,self.z_right, s=20,  cmap="jet", marker="o")

        mp.title("y-z") #图名
        mp.xlabel("y")#x轴标签
        mp.ylabel("z")#y轴标签
        mp.tick_params(axis='both')#x,y轴都有刻度

        mp.show()
    
    def show_rightend(self,xr,yr,zr):
        mp.figure("3D Scatter", facecolor="lightgray")
        ax3d = mp.gca(projection="3d")  # 创建三维坐标

        mp.title('right points', fontsize=20)
        ax3d.set_xlabel('x', fontsize=14)
        ax3d.set_ylabel('y', fontsize=14)
        ax3d.set_zlabel('z', fontsize=14)
        mp.tick_params(labelsize=10)



        ax3d.view_init(elev=0,    # 仰角
                    azim=30       # 方位角
                    )


        ax3d.scatter(self.x_right,self.y_right,self.z_right, s=20,  cmap="jet", marker="o")

        ax3d.scatter(xr,yr,zr, s=20,  cmap='red', marker="o")

        mp.show()

mesh=o3d.io.read_triangle_mesh("cantilever.obj")
obj=find_end(mesh)
obj.show_hull
