# run the PDF generation, the input ntuple is named input_ntuple.root and the generated file is named output.root
# the input_neuple.root file is at the location of (0,0,-500)mm with 2 MeV electron isotropically pointing. Light yield 100 with 1pct wbls.
  ./app2 -i input_ntuple.root\
         -n 100000 \
         -s 0 \
         --sometime 10 \
         --nbins 20 \
         --originZ -500 \
         --ndir 420 \
         --ntime 100 \
         --scanning \
         -r output.root

# once you obtained the pdf, you can run the direction evaluation. 
# There are local minima everywhere so doing a scanning method is more desired.
# fit_ntuple is with 2 mev electron pointing downward at (0,0,-500)mm with ly100 1pct WbLS.
  ./app2 -i fit_ntuple.root\
         -n 10 \  # fit for 10 events
         -s 0 \ # no vertex smearing
         --sometime 4 \ # use a tigher time cut to select more Cherenkov light
         --nbins 20 \
         --originZ -500 \
         --ndir 420 \
         --ntime 100 \
         --scanning \
         --externalPDF pdf_z-500_1pct_ly100.root \
         -o result.txt
