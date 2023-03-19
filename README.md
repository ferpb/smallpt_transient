# smallpm transient

Very simple transient renderer based on Kevin Beason's 99 line path tracer.

Transient rendering breaks the assumption of infinite light speed in conventional renderers, allowing to generate time-resolved images of light propagating through the scene.

The implementation is based on path tracing and uses histogram binning to distribute each sample to the corresponding temporal slice depending on its time of arrival to the sensor.

|![steady](images/steady.png)|![transient delta](images/transient_delta.gif)|![transient continuous](images/transient_continuous.gif)|
|:---:|:---:|:---:|
|Steady state (2e5 samples)|Transient state with delta light pulse (2e5 samples for all frames)|Transient state with continuous lighting (2e5 samples for all frames)|
