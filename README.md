raytracer
=========

CS 775 Ray Tracer in C++ 

***

**Team** 

110050012 - Mayank Meghwanshi (@mayank127)

110050041 - Shivam H Prasad (@shivamh71)

-Wang Xinjiang (swanxinjiang@gmail.com)

***

**Instructions**

To Compile use > Makefile

To Run > ./ray-tracer structure-file

***

**Code Description**

1. **Main File**
 * It contains parser for input scene file.
 * It first takes the light parameters.
 * It then takes camera parameters.
 * Now for each random point from the light a ray is shot ray tracing is done.
 * For ray tracing we find out nearest intersection.
 * Now if specular object then send refracted and reflected rays.
 * If diffused object then light is reemitted from  point of intersection and generate a random emission direction according to (sin^2(theta), phi).

2. **Objects**
 * For each object we find out intersection with a ray from a given origin and direction.
 * We also have a transformation matrix which is used to transform ray to object coordinate system and then use that to find out the intersection.
 * We also have a function which returns normal of the object, which takes as input a point on surface and trasnforms it to object coordinate system and after finding normal at that point transform back to world coordinate space.
 * Different objects made - 
    * **Sphere**
    * **Cylinder** ** Note that the cylinder is put along the z-axis by default **
    * **Cube**

3. **Vec3**
 * This contains function needed for a vector with 3 coordinates x,y,z.
 * It contains functions like dot product, length, cross product, normalize, project, *,/,+,- etc.

4. **Matrix**
 * This contains 3x3 matrix used for rotation in 3-d Space.
 * Also it has usefull functions like inverse, transpose, multiply, transform etc.

***

**Sample Input**

***



**Sources**

* [A sample Ray Tracer for reference](http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-1-writing-a-simple-raytracer/)
