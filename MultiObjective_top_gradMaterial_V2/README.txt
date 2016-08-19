-----------

This is forked from the folder 'SIMP_elastic_heat_gradient'. 

I forked it becasue I want to make a major change to how the design vars are represented and I'm not sure it will reconverge. 

I want 4 or 9 elements to be controlled by single design var. 

The problem is that when I'm using 4 node elements that at the meso scale, I get triangles when trying to maximize the strain energy. I just draws a traingle between the 3 nodes that have displaced the farthest. I'm interested in more complex meso structures, so I need to give a more complex strain field. Hence, I want to have 4 or 9 nodes for each design var. 

Best, yet, If I could make the elements per design var a varriable that could change, then I could have even more control over what is actually going on. 

The goal of this, is to better target how the macor-meso design will look on the macro scale before sending info on the meso scale. If I can better predict what I think the meso scale will generate and the stiffenss macro of the homogenized meso material, then I can go ahead and use this estimated strain energy density on the macro scale. The close I am to the target properties on the first try, then the fewer iteratiosn are needed to stabilize the actual strain energy of each element. 