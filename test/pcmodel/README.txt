# how to

+ install and run the pcdl unit test pcmodel
cp -r ../physicelldataloader/test/pcmodel/ user_projects/
make load PROJ=pcmodel
make
./project

+ set 2d or 3d in config/PhysiCell_settings.xml
<use_2D>true</use_2D>  <!-- 2d -->
<use_2D>false</use_2D>  <!-- 3d -->

+ tar gz the output folders
tar -czvf output_2d.tar.gz output_2d/
tar -czvf output_3d.tar.gz output_3d/
