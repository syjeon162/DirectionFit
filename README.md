A fitting tool for the WbLS detector direction.
There are two kinds of PDF methods : per PMT PDF and universal angle-time PDF.
The input is pre-calculated information in a text file. The reasoning is to remove the ratpac dependency.
An example file is included.

For the formal one, per-calculated PDF will be needed. That calculation can be done with the same tool with the option --perPMT without --externalPDF.
After that's done, both --perPMT and --externalPDF together will activate the reading of external PDFs.

For a lot of features now, you need to ask me why they are there.

An example to run:
$make
$./app2 -i example_fit.txt -m 4 -p pmtlocation_4ton_index_17.txt --evtBase 0 -n 1 -s 0 --sometime 0.2 --upper 1 --lower 0 --orient 999 -b 0 --trueLight --perPMT --nbins 20 --originX 0 --originY 0 --originZ -500 --scanning --externalPDF example_water_perPMT_PDF.root

It will take the extrenal PDF example_water_perPMT_PDF.root, doing per PMT PDF, setting vertex (0,0,-500), taking true light from the file, taking two events.
If you keep --scanning, it will do a scan to find the minimum, otherwise remove it, it will run with Minuit.

```./app2 -i ${INPUT_FILE} -o ${OUTPUT_SCAN} -m ${MASS} -p pmtlocation_4ton_index_${BOTTOM_CASE}.txt -r t ${OUTPUT_ROOT} --evtBase 0 -n 10000 -s 0 --sometime 2 --upper ${UPPER} --lower 0 --orient 999 -b 1 --nbins 20 --originX 0 --originY 0 --originZ ${ZLOCATION} --npmt ${NPMT} --ndir 400 --ntime 30 --perDir --timeCorrection ${TIMECORR}```
Running this will give you pdf
if you specify
--externalPDF ${EXTERNALPDF}
to make it like
./app2 -i ${INPUT} -o ${OUTPUT} -m 4 -p pmtlocation_4ton_index_${BOTTOM_CASE}.txt --evtBase 0 -n 10000 -s 0 --sometime 2 --upper ${upper} --lower ${lower} --orient 999 -b 1 --nbins 20 --originX 0 --originY 0 --originZ ${ZLOCATION} --npmt ${NPMT} --ndir 420 --ntime 30 --perPMT --perDir --externalPDF ${EXTERNAL_PDF} --scanning
then it will load the pre-obtained pdfs. If you don't specify externalPDF, it will only generate pdf and exit, otherwise it will perform the likelihood scanning. 

NPMT, BOTTOM_CASE, MASS, UPPER, ZLOCATION, TIMECORR all depend on the geoemtry
if considering z at (0,0,-500)mm with curved bottom, NPMT=236, BOTTOM_CASE=curved, MASS=4, UPPER=the event number in each of your file, I am taking 5,000 per file (don't use too large, there might be memory issue). ZLOCATION=-500, TIMECORR=0
