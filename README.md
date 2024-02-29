[//]: # (The name of the repository with a brief, high-level description of the
project)
# Pupil Plate Ellipticity (PUPE)- Photometric Analysis Tool (PAT)

## Purpose
The goal of this software tool is to help the user create and analyze pupil
plate images in order to identify the best position for the tip and tilt
of the M2 mechanism during the collimation of a telescope. A "pupil plate"
is a defocused star that is defocused so much that it appears as a donut.
The shape of that donut can be used to characterize (and hopefully be able
to eliminate) optical distortion.

[//]: # 
Users can focus on the following Examples and Usage sections. The rest of the
sections are meant for developers.


[//]: # (Code examples of the usage)
## Examples
Using the engineering password for the "aws site", log onto the
`operational-support` machine via
```
ssh eng@operational-support.aws.lco.gtn
```
The standard offline analysis can be run using the following command:
```
pupe-pat /home/eng/input_data_20180501 /home/eng/output_results_20180501
```
replacing the first path with the directory containing the input images to be
analyzed and the second with the directory to store the output results. Running
```
pupe-pat -h
```
will print the help message with descriptions of the parameters.

The the most up to date docker image and script are automatically installed 
on `operational-support` using Puppet.

To run the real-time analysis, run:
```
pupe-pat-realtime /home/eng/input_data_20180501 /home/eng/output_results_20180501
```
The real-time analysis mode should be run on a machine that has 
/archive mounted, specifically `operational-support`.
Once images are done being taken for the night, press `control-c` to stop
the listener and generate the merged pdfs described below.

[//]: # (A description of how the software in the project is used)
## Usage
During collimation, engineers use a script to take a number of out
of focus images at various tips and tilts of the M2 mechanism. 
Previously, the best image was previously identified by eye.
This code is meant to quantify this process. We currently only support
the LCO 1-m class of telescopes due to scaling the initial guess of the
donut size, but this should be straightforward to extend.

Images to be analyzed should be taken significantly defocused.
The stars should look like a donut with a hole in center. Currently,
we require that the inner hole have a radius of at least 5 pixels. 
The outer edge of the donuts should have a radius of at least 50 pixels
at the moment. For the  Sinistro cameras, this corresponds roughly to a 
defocus of +-12 mm in the focal plane (which is the value stored in FOCDMD 
in the header). There are plans underway to shrink this 
requirement to a 25 pixel radius which would allow defocus parameters down 
to +-6 mm. This is currently a requirement in pixels rather than on defocus,
so for the SBIGs in 1x1 binning we can currently get down to +-7.5 mm. 

The pupil plate images need to be taken on a bright, isolated star. A cut of 
V~4 seems to work well. This produces images that are high signal-to-noise
(half the saturation value) in about 1 second exposures. Fainter stars can be
used if they are isolated (the donuts do not overlap) and the exposure times
are long enough to get ~20,000 counts in the images.

The code automatically detects sources using Source Extractor in Python,
SEP, https://sep.readthedocs.io/en/v1.0.x. The code then makes a cutout
around the donut and fits a model using astropy.modeling that has one 
value in the "hole" of the donut, one brightness for the pupil plate and 
 one for the background level. The centers and radii of the circles 
that form the boundaries of the pupil plate are also allowed to vary. 
There is code in this repo to fit ellipses instead of circles and to fit
brightness gradient across the donut, but they are currently disabled.
Furthermore, non-circular sources are filtered based upon ellipticity if
the A/B-axis ratio (as reported by `sep.extract`) exceeds 1.15.

The results are saved in the output directory specified by the user at runtime.
Each star in each image that is fit as a pupil plate produces a pdf file
with the filename of the image followed by a source id. At the end of the
analysis the pdfs are merged by defocus position into
pupe-pat-{focdmd}-donut.pdf (default) where {focdmd} is replaced 
with the requested focus of the images. The code also produces a "quiver"
plot (replacing donut with quiver in the name above). This shows the magnitude
and direction of the offset of the centers of the inner and outer edges of
the pupil plate. The "best" M2 tip and tilt are marked with a red X.

The code can be run in two modes (see below): real-time and offline modes.
The offline mode will analyze all images in a given directory
that are EXPOSE or EXPERIMENTAL obstypes, and are defocused (FOCDMD != 0.0).
By default, `PROPOSAL_ID` is not used to filter images. However, the
`--proposal-id PROPOSAL_ID` argument may be specified on the command line
to analyze only those images that are taken by the specified `PROPOSAL_ID`.
To ensure that the proposal is correctly included in the header, the Sequencer
should be set to manual when taking the images.

The real-time mode has the same criteria for which images to analyze, but
watches the input directory for new files. Right before an M2 tip/tilt run,
this watcher should be started on the directory where the files will appear. 
For the SBIG cameras, this can be the `raw` directory for the given DAY-OBS.
The Sinistro frames need the amplifiers to be stitched together before they
can be analyzed, so the watcher should be directed to the `preview` directory
instead.

[//]: #
## Installation
The code is currently installed via Puppet. The only machine that has the code
automatically installed is `operational-support`.

On other machines, the code is meant to be run in a docker file.
You can pull the docker image by typing the following: 
```
docker pull docker.lco.global/pupe-pat:{0.1.0}
```
replacing {0.1.0} with the current version. You can find the current version 
checking <https://github.com/lcogt/pupe-pat/releases>.

Alternatively, you can install the code directly 
(after installing the dependencies, hopefully in an environment)
by running
```
git clone https://github.com/lcogt/pupe-pat
cd pupe-pat
python setup.py install
```
but the docker installation is preferred. This setup will not include the wrapper
scripts described above in the Examples section above. For using a 
machine besides `operational-support` or using a direct python installation
see the Advanced Examples below.

### Prerequisites
Docker needs to be installed on your system and you need to be a docker user.
`eng` is probably the best option. `operational-support` is the only machine that
will have the pupe-pat and pupe-pat-realtime wrapper scripts installed via 
Puppet currently.

### Dependencies
- sphinx

- numpy

- astropy

- lcogt_logging

- sep

- matplotlib

- watchdog

- scipy

- emcee (only needed for manual analysis)

- corner (only needed if emcee is installed)

- pypdf2


[//]: # (Describe how to run tests in the project)
## Tests
Currently, the tests are run manually. Our test data set was taken on
kb05 at bpl on 20180409. We use the whole raw directory of images.
To reduce the test data run the following on the `operational-support` machine
as the `eng` user:
```
pupe-pat /home/eng/pupe-pat-test-data /home/eng/pupe-pat-tests-results
```
You can test other versions using the docker commands
in the Advanced Examples section below.

The tests pass if all the blue circles align with the donuts
in the PDF files produced in the output directory.

Eventually this procedure will be replaced with automatic testing in Jenkins.

[//]: # (Details on how to deploy the project; if the project is deployed using
CD, just put a link to the job in this section)
## Deployment
Currently, new versions are built manually. You can run
```
docker build -t docker.lco.global/pupe-pat:{0.1.0} .
docker push docker.lco.global/pupe-pat:{0.1.0}
```
replacing {0.1.0} with the current version. You can find the current version 
checking <https://github.com/lcogt/pupe-pat/releases>.
Eventually we will replace this procedure with continuous deployment in Jenkins.

After the updated version of the docker image has been pushed, update and
deploy the puppet config files to use the new version.

[//]: #
## Advanced Examples
On machines other than `operational-support`, if the script is installed in a python
environment, the standard offline analysis can be
run using the following command:
```
run_pupepat --input-dir /home/eng/input_data_20180501 --output-dir  /home/eng/output_results_20180501
```
replacing the paths with the desired input and output paths.

The other available parameters are `--output-table` and `--output-plot`, and `--proposal-id`.
These control the names of the output files produced at the end of the analysis.
For example, if you wanted to run the analysis a second time but not overwrite
the results, you could run
```
run_pupepat /home/eng/input_data_20180501 --output-dir  /home/eng/output_results_20180501 --output-table pupe-pat-run2.dat --output-plot pupe-pat-run2
```
Note that the `--output-table` flag needs an extension (e.g. `.dat` or `.txt`),
but the `--output-plot` flag does not. It automatically makes pdf files. 
This difference is because each focus position produces a separate plot.
The `--proposal-id` flag can be used to set which proposal took the pupil plate images. 
The default is to not filter on proposal id and process all images.

The real-time analysis has identical arguments, but is run as
```
run_pupepat_listener /home/eng/input_data_20180501 --output-dir  /home/eng/output_results_20180501
```

On machines other than `operational-support` the preferred method to run the code
is via docker. You can do this by running
```
docker pull docker.lco.global/pupe-pat:{0.1.0}
docker run --rm -v /home/eng/pupe-pat-test-data:/input -v /home/eng/pupe-pat-tests-{0.1.0}:/output docker.lco.global/pupe-pat:{0.1.0} run_pupepat_listener --output-dir /output --input-dir /input
```
replacing {0.1.0} with version of interest. You can find the current version 
by checking <https://github.com/lcogt/pupe-pat/releases>. The mount points inside the container are arbitrary names, needing only to match
the `--input-dir` and `--output-dir` arguments.

Running the real-time analysis using docker is similar
```
docker pull docker.lco.global/pupe-pat:{0.1.0}
docker run --rm -it -v /home/eng/input_data_20180501:/input -v /home/eng/pupe-pat/output_results_20180501:/output docker.lco.global/pupe-pat:{0.1.0} run_pupepat_listener --output-dir /output --input-dir /input
```
The real-time analysis requires the -it flags so that it will receive the 
`control-c` command to stop the listener and make the final set of pdfs.

[//]: #
## Acronyms

| Acronym | Description                                |
|---------|--------------------------------------------|
| M2      | Secondary Mirror mechanism on 1m telescope |
| FWHM    | Full width half max                        |

[//]: #
## License
This code is licensed under GPL v3.0. See LICENSE for more information.

[//]: #
## Support
 [Create an issue](https://issues.lco.global/)

---------------------------------------------------------------
