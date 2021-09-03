import numpy as np
import open3d as o3d
import plotly.graph_objects as go


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
        import plotly.graph_objects as go


    """
    
    
    def __init__(self,mesh):
        self.pcd=mesh
        self.point=np.asarray(self.pcd.vertices)

    def get_hull(self):
        self.hull, _ = self.pcd.compute_convex_hull()
        self.hull_ls = o3d.geometry.LineSet.create_from_triangle_mesh(self.hull)
        self.point2=np.asarray(self.hull.vertices)
    
        self.hull_ls.paint_uniform_color((1, 0, 0))
        o3d.visualization.draw_geometries([self.pcd, self.hull_ls])
    
    def show_points(self):
        x=self.point2[:,0]
        y=self.point2[:,1]
        z=self.point2[:,2]

        fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z,
                                   mode='markers',
                                   marker=dict(
                                    size=5,
                                                    # set color to an array/list of desired values
                                    colorscale='Viridis',   # choose a colorscale
                                    opacity=0.8
                                     ))])

        fig.show()


mesh=o3d.io.read_triangle_mesh("cantilever.obj")
obj=find_end(mesh)
obj.get_hull()
obj.show_points()