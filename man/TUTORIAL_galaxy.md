# PhysiCell Data Loader Tutorial: pcdl and Galaxy

[Galaxy](https://en.wikipedia.org/wiki/Galaxy_(computational_biology)) is an open-source, web-based platform for bioinformatics analysis.

There a hand full of galaxy server instances around the world (USA, Europe, France, Australia, possibly others).
We here focus on [the European instance of Galaxy](https://usegalaxy.eu/) since it has pcdl tools installed.

Special thanks to Björn Grüning, who helped me to port pcdl into Galaxy!


## &#x2728; Run PhysiCell Studio on Galaxy and save the output

1. Point your web browser to the European instance of Galaxy: https://usegalaxy.eu/
1. Required! Log into your user account. Create one if you do not already have one.
1. On the left-side panel, click "Interactive Tools".
1. Near the left-side top corner, in the "Tools" search bar, type "physicell studio".
1. Choose "PhysiCell Studio".
1. On the center panel, click on one of the "Run Tool" buttons.
1. On the left-side panel, after a while, "Active Interactive Tools" should pop up.
1. Under "Active Interactive Tools" click "PhysiCell Studio". 
1. Run a simmulation: In brief, in the middle panel, click the "Run" tab, then click the "Run simmulation", then wait until on the screen output "Processed finish" shows up.
1. On the center panel menu bar, choose "Misc > put on History > all output.zip".
1. On the right-side History panel, after a while, an "all\_output.zip" file should pop up.

### Unzip the zip folder
1. Near the top left corner, click "Tools".
1. Near the left-side top corner, in the "Tools" search bar, type "unzip".
1. Choose "Unzip Unzip a file".
1. Mouse drag and drop the "all\_output.zip" file from the right-side History panel into the tool's parameters "Input file" field in the middle panel.
1. On the center panel, click on one of the "Run Tool" buttons.
1. On the right-side History panel, after a while, an unzipped file collection should pop up.


## &#x2728; Upload local PhysiCell output to Galaxy

1. Point your web browser to the European instance of Galaxy: https://usegalaxy.eu/
1. Optionally: Log into your user account.
1. Near the top left corner, click "Upload (New upload Beta)".
1. Near the top left corner, choose "Upload from Computer".
1. In the center panel click "Brows Files".
1. Click your self through to the PhysiCell output folder.
1. Click on one file in the output folder, then, on your keyboard, press Ctrl + A to choose all files.
1. Click "Open". After a while the filenames will appear in the central panel.
1. On the center panel bottom, flip the switch "Create a collection from these files" on.
1. "Enter a collection name" (e.g. output) and choose  "Collection Type" "List".
1. Click "Start".
1. After a while, the uploaded data collection should appear in the in the right-side History panel.


## &#x2728; The pcdl tools on Galaxy

1. Point your web browser to the European instance of Galaxy: https://usegalaxy.eu/
1. Optionally: Log into your user account.
1. Upload or generate PhysiCell output as described above.
1. Near the top left corner, click "Tools".
1. Near the left-side top corner, in the "Tools" search bar, type "pcdl".
1. Choose one of the pcdl\_ tools (e.g. pcdl\_get\_version).
1. Mouse drag and drop the uploaded or generated PhysiCell output data collection in the left-side History panel into the "essential:"  "data collection" field in the middle panel.
1. Optionally, tweak the "essential" and "advanced" parameters.
1. Click "Run Tool".
1. After a while, the result file(s) should appear in the left-side History panel.


## &#x2728; Further readings

You can learn more about Galaxy here:
+ https://galaxyproject.org/tutorials/g101/

How to load, develop, or run a PhysiCell model in the studio can be learned here:
+ https://physicell-studio.readthedocs.io/en/latest/index.html
+ https://github.com/PhysiCell-Tools/PhysiCell-Studio

That's it. The rest is analysis within Galaxy!
