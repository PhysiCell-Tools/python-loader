# PhysiCell Data Loader Tutorial: pcdl and Galaxy

[Galaxy](https://en.wikipedia.org/wiki/Galaxy_(computational_biology)) is an open source, web-based platform for bioinformatics anlysis.

There a hand full of galaxy server instances around the world (USA, Europe, France, Australia, possibly others).
We here focus on [the European instance of Galaxy](https://usegalaxy.eu/) since it has pcdl tools installed.

Special thanks to Björn Grüning, who helped me to port pcdl into Galaxy!


## &#x2728; Run PhysiCell Studio on Galaxy

1. In your web browser, open the European instance of Galaxy: https://usegalaxy.eu/
1. Required, log into your user account. Create one if you do not already have one.
1. One the left side pannel, click "Interactive Tools".
1. Near the left-side top corner, in the "Tools" search bar, type "PhysiCell Studio".
1. Choose "PhysiCell Studio".
1. On the center pannel under additional option flip the switch "Email notification" to Yes. 
1. On the center pannel click on one of the "Run Tool" buttons.
1. On the left-side pannel, after a while, "Active Interactive Tools" should pop up.
1. Under "Active Interactive Tools" click "PhysiCell Studio". 
1. How to load, develop, or run a PhysiCell model in the studio can be learned here:
   + https://physicell-studio.readthedocs.io/en/latest/index.html
   + https://github.com/PhysiCell-Tools/PhysiCell-Studio
1. In the left side "History" pannel do NOT DELETE the "PhysiCell Studio on : all reults archive."
1. Run simmulation.
1. [TKBue]



## &#x2728; Upload local PhysiCell output to Galaxy

1. In your web browser, open the European instance of Galaxy: https://usegalaxy.eu/
1. Optionally, log into your user account. Create one if you do not already have one.
1. Near the top left corner, click "Upload" (New upload Beta).
1. Near the top left corner, choose "Upload from Computer".
1. In the center panel click "Brows Files".
1. Click your self through to the PhysiCell output folder.
1. Click on one file in the output folder, then on your keyboard press Ctrl + A to choose all files.
1. Click "Open". After a while the filenames will appear in the central panel.
1. On the center panel bottom, flip the switch "Create a collection from these files" on.
1. "Enter a collection name" (e.g. output) and choose  "Collection Type" "List".
1. Click "Start".
1. After a while, the uploaded data collection should appear in the history in the right-side panel.


## &#x2728; The pcdl tools on Galaxy

1. In your web browser, open the European instance of Galaxy: https://usegalaxy.eu/
1. Optionally, log into your user account. Create one if you do not already have one.
1. Upload or generate PhysiCell output as described above.
1. Near the top left corner, click "Tools".
1. Near the left-side top corner, in the "Tools" search bar, type "pcdl".
1. Choose one of the pcdl_ tools (e.g. pcdl\_get\_version).
1. Mouse drag and drop the uploaded or generated PhysiCell output data collection in the left-side history panel into the "essential:"  "data collection" field in the middle panel.
1. Optionally, tweak the "essential" and "advanced:" parameters.
1. Click "Run Tool".
1. After a while, the result file(s) should appear in the left-side history panel.


## &#x2728; Further readings

You can learn more about Galaxy here:
+ https://galaxyproject.org/tutorials/g101/

That's it. The rest is analysis within Galaxy!
