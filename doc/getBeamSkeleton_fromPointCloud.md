
<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [Package: trimesh](#package-trimesh)
  - [1. Load mesh from  .stl](#1-load-mesh-from--stl)
  - [2. Functions of mesh](#2-functions-of-mesh)
  - [3. Get cross section](#3-get-cross-section)
- [Class: SkeletonOfBeam](#class-skeletonofbeam)
  - [1. Using process](#1-using-process)
  - [2. Source code](#2-source-code)

<!-- /code_chunk_output -->




# Package: trimesh

```py
import trimesh
```

## 1. Load mesh from  .stl
.obj file 在导入的时候会有些问题，此时需要先转成 .stl
```py
mesh = trimesh.load_mesh("cantilever.stl")
```

## 2. Functions of mesh
Centroid
```py
[1*3 array] mesh.centroid
```

XYZ Bounds:
```py
[2*3 array] mesh.bounds

>>>
[[xmin, ymin, zmin]
 [xmax, ymax, zmax]]
```

## 3. Get cross section

> **3D slice to 2D slice**
```py
# 3D slice
slice_3D = mesh.section(plane_origin=point, plane_normal=vector)
# 2D slice
slice_2D, to_3D = slice_3D.to_planar()
```

> **Properties of slices**
```py
[1*3 array] slice_2D.centroid
[double]    slice_2D.area
```

> **View**
```py
slice_3D.show()
slice_2D.show()
```

> **Split slice (划分截面数据，可用于初步降噪)**
```py
[trimesh.Path2D[]]
slices_splited = slice_2D.split()
```
这时原先的 `slice_2D` 已经被划分为几个小的子面（可以理解为 cluster），因此可以通过选取面积最大的子面，实现初步降噪
```py
sliceIndex = np.argmax([s.area for s in slices_splited])
slice_2D = slices_splited[sliceIndex]
```
大致效果如下
<center>
    <img src="https://github.com/ZhekaiLi/PICTURE-for-markdown/raw/master/2021-07/pic_202107301655.jpg" style="zoom:80%"> <img src="https://raw.githubusercontent.com/ZhekaiLi/PICTURE-for-markdown/master/2021-07/pic_202107301656.jpg" style="zoom:80%"> <br>
    <div style="color: #999;"></div>
</center><br>

> **旋转摆正**
```py
[3*3 array] slice_2D.apply_obb()
slice_2D.show()
```
<center>
    <img src="https://github.com/ZhekaiLi/PICTURE-for-markdown/raw/master/2021-07/pic_202107301657.jpg" style="zoom:0%"> <br>
    <div style="color: #999;"></div>
</center><br>

> **Several cross sections**
```py
sections = mesh.section_multiplane(plane_origin=point, 
                                    plane_normal=vector,
                                    heights=heightsList)
sections = [i for i in sections if i is not None]
```
对每个 section 都是一个 slice_2D，同样可以对其进行各种基础操作、数据读取

还可以可以将它们合并后展示，总体流程：
```py
z_extents = mesh.bounds[:,0]
z_levels  = np.arange(*z_extents, step=2) # 每隔两个单位长度截一次
sections = mesh.section_multiplane(plane_origin=mesh.centroid, 
                                    plane_normal=[1, 0, 0],
                                    heights=z_levels)
sections = [i for i in sections if i is not None]
combined = np.sum(sections)
combined.show()
```
<center>
    <img src="https://github.com/ZhekaiLi/PICTURE-for-markdown/raw/master/2021-07/pic_202007301714.jpg" style="zoom:0%"> <br>
    <div style="color: #999;"></div>
</center><br>


-----------------------


# Class: SkeletonOfBeam
## 1. Using process
```py
# Create skeletonOfBeam object from a trimesh.mesh and a rough skeleton vector
sob = SkeletonOfBeam(mesh, rough_normalVector)

# Get intersections first, then get centroids of intersections
sob.getIntersections(step=1)
sob.getCentroids()

# Read the properties
[np.array] sob.centroid
[trimesh.Path3D[]] sob.Intersections
[np.array[]] sob.Centroids
```


## 2. Source code 
```py
class SkeletonOfBeam:
    """
    Get skeleton of beam

    Arributes:
        mesh: Input trimesh.mesh
        centroid: An array of the coordinates of mesh's centroid
        
        Intersections: A list of trimesh.Path3D
        Centroids: A list of np.array
    """
    mesh = None
    centroid = None
    Intersections = None
    Centroids = None
    
    def __init__(self, mesh, rough_normalVector):
        self.mesh = mesh
        self.centroid = mesh.centroid.copy()
        # 质心截面的一个大致的法向量（目前单指 [1, 0, 0]）
        self.nVec = rough_normalVector
        
          
    def getIntersections(self, step=1):
        """Get the intersections of beam along vector [1, 0, 0]
        
        :param step: interval of intersections
        """
        sections = []
        extents = mesh.bounds[:, 0] # 截取区间
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
                    sections.append(slice_2D.to_3D(to_3D))
            except:
                pass
        
        self.Intersections = sections
                
    def getCentroids(self):   
        """Get the centroids of intersections
        """
        self.Centroids = []
        for s in self.Intersections:
            self.Centroids.append(s.centroid)    
    
            
    def returnXYZOfCentroids(self):
        return np.array(self.Centroids)[:, 0], np.array(self.Centroids)[:, 1], np.array(self.Centroids)[:, 2]
```

