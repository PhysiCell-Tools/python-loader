# PhysiCell Data Loader Tutorial: pcdl and Galaxy

[Galaxy](https://en.wikipedia.org/wiki/Galaxy_(computational_biology)) is an open source, web-based platform for bioinformatics anlysis.

There a hand full of galaxy server instances around the world (USA, Europe, France, Australia, possibly others).
We here focus on [the European instance of Galaxy](https://usegalaxy.eu/) since it has pcdl tools installed.

Special thanks to Björn Grüning, who helped me to port pcdl into Galaxy!


## Run PhysiCell Studio on Galaxy

1. In your web browser, open the European instance of Galaxy: https://usegalaxy.eu/
2. Optionally, log into your user account. Create one if you do not already have one.
3. Near the top left corner, click "Tools".
4. Near the left-side top corner, in the "Tools" search bar, type "PhysiCell Studio".
5. Choose "PhysiCell Studio".
6. Choose "Run Tool" (center panel or top right).
7. How to load, develop, or run a PhysiCell model in the studio can be learned here:
   + https://physicell-studio.readthedocs.io/en/latest/index.html
   + https://github.com/PhysiCell-Tools/PhysiCell-Studio
8. [TKBue]


## Upload local PhysiCell output to Galaxy

1. In your web browser, open the European instance of Galaxy: https://usegalaxy.eu/
2. Optionally, log into your user account. Create one if you do not already have one.
3. Near the top left corner, click "Upload" (New upload Beta).
4. Near the top left corner, choose "Upload from Computer".
5. In the center panel click "Brows Files".
6. Click your self through to the PhysiCell output folder.
7. Click on one file in the output folder, then on your keyboard press Ctrl + A to choose all files.
8. Click "Open". After a while the filenames will appear in the central panel.
9. On the center panel bottom, flip the switch "Create a collection from these files" on.
10. "Enter a collection name" (e.g. output) and choose  "Collection Type" "List".
11. Click "Start".
12. After a while, the uploaded data collection should appear in the history in the right-side panel.


## The pcdl tools on Galaxy

1. In your web browser, open the European instance of Galaxy: https://usegalaxy.eu/
2. Optionally, log into your user account. Create one if you do not already have one.
3. Upload or generate PhysiCell output as described above.
4. Near the top left corner, click "Tools".
5. Near the left-side top corner, in the "Tools" serach bar type "pcdl".
6. Chosse one of the pcdl_ tools (e.g. pcdl\_get\_version).
7. Mouse drag and drop the uploaded or generated PhysiCell output data collecton in the left-side history panel into the "essential:"  "data collection" field in the middle panel.
8. Optionally, tweak the "essential" and "advanced:" parameters.
9. Click "Run Tool".
10. After a while, the result file(s) should appear in the left-side history panel.


### &#x2728; Further readings

You can learn more about Galaxy here:
+ https://galaxyproject.org/tutorials/g101/

That's it. The rest is analysis within Galaxy!
