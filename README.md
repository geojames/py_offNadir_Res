# py_offNadir_Res

**Name:**          offNadir_resolution_v1-0.py

**Compatibility:**  Python 3.5

### Description:
This program calculates the pixel resolution and instanteous
field of view for a camera from an input file provided by
the user. The file nees to be a comma-delimited text file 
(*.CSV) with a header row:
                    
                    Name,focal,sensor_x,sensor_y,pixel_x,pixel_y,flyH,angle
                    P3,3.61,6.24,4.71,6000,4000,100,10
                    
where...
- **Name** is camera name
- **focal** is the focal length of camera (mm)
- **sensor_x** is the sensor's long dimension (mm)
- **sensor_y** is the sensor's short dimension (mm)
- **pixel_x** is the number of pixels along the long dimension
- **pixel_y** is the number of pixels along the short dimension
- **flyH** is the flying height/altitude of the camera
- **angle** is the off-nadir angle of the camera (0 = nadir, 90 = horizontal)


*the software will only work for low-oblique images (i.e. no horizon visable) and will spit out a error if the off-nadir angle puts the top part of the verical field of view over the horizon.

### To Run:
- Run the software from a Python shell or editor
	- A file chooser will pop-up (it may try to hide)
- Choose your CSV file
- The program will run and create a formated text file with the output in the same location as the input file. The name will be the same with the suffix '_resolution'
