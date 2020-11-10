This folder is used to determine paleo-locations of sediment sites.

Execute the following steps:

1. Put the sites of which you would like to determine the paleo-location in the dictionary 'sites' in 'write_shapefile.py' .

2. Run 'write_shapefile.py' , then the shapefile 'sitelocs' is written. Use the shapefile.yml conda environment (or installing the 'pyshp' package should probably be enough). 

3. Open the rotation file 'vanHinsbergen_master.rot' and the 'sitelocs' files in GPlates (and if you like the Coastlines.gpml for orientation;'File'->'Manage Feature Selections'->'Open File').

4. Go from 0Ma to the time period of interest (38Ma in this case). Export a shapefile from Gplates which contains the rotated paleo-locations of the sediment sample sites. 

5. Write the sitelocs as shapefile in Gplates:
	- Make sure all of the outlines are set to visible i.e. they have a little eye checked
	  in the Layers dialog.
	- Choose 'Reconstruction', 'Export...' from the general GPlates menu (top of screen);
	- Select 'Export Single Snapshot Instant' and click 'Use Main Window Type', or just
	  enter '38' into the 'Time' box.
	- Change 'Target Directory' into the 'Geometries' subfolder within 'Adjusted Files'.
	- Click 'Add Export', 'Reconstructed Geometries' and choose 'Shapefiles (*.shp)'.
	- Make sure that 'Export to a single file' is unchecked and 'Export to multiple files'
	  (as well as 'Separate output directory') is checked. 'Wrap polyline and polygon
	  geometries to the dateline' should also be checked to avoid horizontal lines.
	- Leave the export name as it is (since it will be stored in separate folders), apart
	  from changing '_%0.2fMa' to '_%0.0fMa'. This will give you the suffix '38Ma' instead
	  of '38.00Ma', which is unnecessary.
	- Click OK and wait for the export to finish (should take just seconds).

6. The new shapefile contains the paleolocations. 'readshapefile.py' is a python script which can be used to read this shape file, and writes the coordinates that can be used further in the analysis.





