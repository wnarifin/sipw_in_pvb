Steps:

In terminal,

1. Make directories

    > `mkdir exp1_2  exp2_2  exp3_2  exp4_2  exp5_2  exp6_2`

1. Copy run.R to each folder

    > `xargs -n 1 cp -v run.R<<<"exp1_2  exp2_2  exp3_2  exp4_2  exp5_2  exp6_2"`

2. Copy par.R to each folder

    > `xargs -n 1 cp -v par.R<<<"exp1_2  exp2_2  exp3_2  exp4_2  exp5_2  exp6_2"`

    then edit according to exp_all_settings.txt

3. Run runall.sh
