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
