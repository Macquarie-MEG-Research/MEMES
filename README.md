# MEMES
**MRI Estimation for MEG Sourcespace (MEMES)** is a set of tools for estimating an appropriate structural MRI for MEG source analysis in Fieldtrip and/or SPM. 

These set of scripts are customised for data acquired from the Macquarie/KIT MEG laboratory using a 160-channel Yokogawa MEG system.

![MEMES](https://github.com/Macquarie-MEG-Research/MEMES/blob/master/actual_memes/3memes.png)

### Outline

MEMES is based on the approach of [Gohel et al., (2017)](https://www.frontiersin.org/articles/10.3389/fninf.2017.00050/full). It uses an Iterative Closest Point (ICP) algorithm to match participant's headshape information to a database of 95 template MRIs from the [Human Connectome Project (HCP)](https://db.humanconnectome.org/app/template/Login.vm;jsessionid=328FAED54C595F8EB6897003A5526C84)

The code then selects the most appropriate MRI and creates a 3D source-model from this for subsequent source analysis.

### Results

So far the results seem promising, with only very minor differences to the patterns/peaks of source localistion compared to real anatomical MRIs

![Results](https://github.com/Macquarie-MEG-Research/MEMES/blob/master/actual_memes/results_1.png)

### Things to work on:

- Validate with larger dataset
- Ground truth?
- Average over first 10 (or so) best fitting MRIs?

