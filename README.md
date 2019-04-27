# State-dependent effects of electrical stimulation on populations of excitatory and inhibitory neurons

Electrical stimulation of neural populations is a key tool for understanding neural dynamics and developing treatments. To investigate the effects of external stimulation on a basic cortical motif, we analyse the dynamical properties of an efficient mean-field neural mass model of excitatory and inhibitory adaptive exponential integrate-and-fire (AdEx) neurons and validate the results using detailed network simulations. 
The state space of the mean-field model and of the detailed network are closely related. They reveal asynchronous up and down-states, bistable regions, and oscillatory regions corresponding to fast excitation-inhibition and slow excitation-adaptation feedback loops. 
Within this dynamical landscape, external stimuli can cause state transitions, such as turning on and off oscillations. Oscillatory input can frequency-entrain and phase-lock endogenous oscillations. The effects of external stimulation are well-predicted by the mean-field model, further underpinning the utility of low-dimensional neural mass models. 	

## Getting Started

### Prerequisites

Please install [anaconda](https://www.anaconda.com/distribution/) for python 2.7 on your computer before you clone this repository.

### Installing

To make it easier to run the code in this project and reproduce our findings, we have created an anaconda environment. This helps you to get all the correct versions of the libraries we have used as well as python 2.7. 

The following command creates the anaconda environment on your computer:

```
conda env create -f nip27.yml
```

## Tutorial

[tutorial]

## Built With

* [Brian2](https://github.com/brian-team/brian2) - A clock-driven simulator for spiking neural networks

## Authors

Caglar Cakan*, Klaus Obermayer
Department of Software Engineering and Theoretical Computer Science, Technische Universität Berlin, Germany
Bernstein Center for Computational Neuroscience Berlin, Germany

*cakan@ni.tu-berlin.de

## License

This project is licensed under the BSD 2 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

* A description of the mean-field population model used here can be found [here](https://github.com/neuromethods/fokker-planck-based-spike-rate-models) as described in: Augustin*, Ladenbauer*, Baumann, Obermayer, Low-dimensional spike rate models derived from networks of adaptive integrate-and-fire neurons: comparison and implementation, PLOS Computational Biology 2017
* The code for precomputing the filter components of the linear-nonlinear cascade is based on this project as well.
