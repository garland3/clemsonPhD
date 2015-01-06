
union(){
 import("C:\\Users\\Anthony G\\SkyDrive\\GarlandINC\\3DPrintingAsAService\\Orders\\ShariBox\\Part2.STL");
 
 

translate([30, 20, 25]){
    linear_extrude(height = 4, center = true, convexity = 2, twist = 0)
   text("Shari", font = "Liberation Sans",halign ="center");
}
 }
 
