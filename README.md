# ChIPAnalysis

Make sure to do

chmod a+x ProcessBAMFile.R

Now you should be able to run

[PATH]/ProcessBAMFile.R --help

to get the detailed usage.

To run basic processing on a BAM file and save to data

[PATH]/ProcessBAMFile.R -m -c -q -n -T --bamdir=BAM --datadir=Data --metadir=Meta --QCdir=QC --trackdir=Tracks --maxfraglen=220 --minfraglen=50 --centerwidth=50 [NAMES OF BAM FILES]

This assumes that your BAM files are in directory "BAM", and that you
have directories "Data", "QC", "Meta", and "Tracks" for saving
different file types. This commenad will do a lot of the work and can
be run in parallel for smaller groups of files.

If you want to generate QC summary then

[PATH]/ProcessBAMFile.R --sizes=50,220,1000 -S QCsummaryFile.csv

the "--sizes" part specifices the fragment sizes you want to focus on.

To generate nuclesomes table (as RData)

[PATH]/ProcessBAMFile.R -N NucSummary.rdata

To generate nuclesomes sheet (as CSV)

[PATH]/ProcessBAMFile.R --Nucsummary NucSummaryFile.csv


