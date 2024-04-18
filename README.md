A fitting tool for the WbLS detector direction.

There are three kinds of PDF methods : per direction PDF, per PMT PDF and universal angle-time PDF.  
Unless you have a very strong reason, we should use the per direction PDF since that is the most precise method to evalute the direction fitting performance by far.  
A bit more information about the fitter itself: https://nino.lbl.gov/eos/DocDB/cgi/ShowDocument?docid=139  

The fitting itself can be done in two ways:  
- minuit minimizer
- scanning (this way is preferred. There will be a lot of local minimer if you do minuit.)

The vertex fitting can be performed simultaneously but the framework needs to be updated to include that feature. At the moment, fix the vertex location.  
The simultaneous fit will take more memory and reduce the running speed so for fast evaluation of the fixed-position sources such as beta, gamma, use this tool is much more efficient. Fast results should be provided to 1. optimize the source design 2. fine tune the framework itself.  

Three useful files are provided to start an example run:
- input_ntuple.root: isotropic 2 mev electrons at (0,0,-500) mm in ly100 1pct wbls.
- fit_ntuple.root: downward pointing 2 mev electrons at (0,0,-500) mm in ly100 1pct wbls.
- pdf_z-500_1pct_ly100.root: pre-calculated PDF with 2,000,000 events generated with isotropic 2 mev electrons at (0,0,-500) mm in ly100 1pct wbls.

There are two steps to complete a fitting procudure, the pdf generation and the direction evaluation.  
A fast way to start is to look at `testRun.sh` and try to run that file. Following provides a bit more information:

### PDF generation  
The common useful information necessary:
- input ntuple file
- number of events
- vertex smearing in mm
- number of bins in theta angle (use 20 if no strong preference)
- original z location in mm
- number of directions (should be larger than number of beta ^ 2, for nbins for theta = 20, use at least 400 for this)
- number of bins in the time grid. Each step is 0.2 ns (use 100 if no other strong preference).
- output file

An example of running this is:
```
./app2_dic -i input_ntuple.root\
     -n 100000 \
     -s 0 \
     --sometime 10 \
     --nbins 20 \
     --originZ -500 \
     --ndir 420 \
     --ntime 100 \
     --scanning \
     -r output.root
```
         
### Direction evaluation  
The common useful information is similar the PDF generation. However, to distinguish the pdf generation and actual fitting. Following things are important:
- name an external pdf
- have a output text file ready

An example of running the result:
```
./app2_dic -i fit_ntuple.root\
     -n 10 \  # fit for 10 events
     -s 0 \ # no vertex smearing
     --sometime 4 \ # use a tigher time cut to select more Cherenkov light
     --nbins 20 \
     --originZ -500 \
     --ndir 420 \
     --ntime 100 \
     --scanning \ # using the scanning method
     --externalPDF pdf_z-500_1pct_ly100.root \ # need to provide per-generated pdf
     -o result.txt # output file
```

### Output
The output file format will be (a,b,c,d):
- a: event id
- b: theta angle (zenith) in the unit of angle = $b * (\pi/N_{bin})$. The upward pointing is 0 and the downward pointing is 1. If you have a downward event, you might expect the b at 19 is better.
- c: similar to the zenith angle, this is the azimuth angle. angle = $c * (2\pi/N_{bin})$.
- d: likelihood value. Smaller value indicates a more likely location.

if you see something like (a,b,-1,c), this only appear once for each event at the end of the event list. It tells you the best fit directly:
- a: zenith
- b: azimuth
- c: marker -1
- d: likelihood

By finding the mark -1, you can plot the overall reconstructed direction.  
Once again, the z vertex can be fitted by yourself with a list of different z. Once you combine all the result list and try to find the minimum value across both direction and z, you can effectively fit the z simultaneously. Imbedding z feature significantly slow the pdf producion so it is in the process of optimizing this procedure.  

### Updates
Update on Jun 13 2023:  
> I have updated the capability of dichroicon and directional source container.
The app2 excutable may not be compatible. Change to app2_dic instead.  
If want to use dichroicon, add the command `--dic`;  
If want to use directional source container, add the command `--source`.

Update on Apr 17 2024:  
> When using command `--source`, the source direction is now retrieved from
the metadata from the eosntuple file.  
Set `/rat/proclast eosntuple` to get an eosntuple output from EosSims.  
