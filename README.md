
# EEG Toolkit

> :last_quarter_moon_with_face: Developed at the [Laboratory for the Study of
> Sleep Slow-wave activity](https://www.med.upenn.edu/slowwavelab/)

A scientific package for computational EEG analysis with an emphasis on 
methodological transparency. Current features:

- EEG Computational Toolkit
    - Loading EEG data
    - EEG visualization
    - Sleep stage handling
    - Power spectral analysis

# Example

I have a test EDF file with a full-night EEG called `edf2.edf`.

```julia
eeg = EEG("edf2.edf", 500, 30, [])
signal = epoch(eeg, 100, 200, "EEG C3-A2")
S = Spectrogram(signal, 500, 3, 0.5)
p = plot_spectrogram(S, 30.0, 2)
```

![Image]("imgs/spetrogram_plot.png")

We can also take a look at the spectrogram as a heatmap:

```julia
p = plot_spectrogram(S, 30.0, 1, :inferno) # Color scheme inferno is better for heatmaps
```

![Image]("imgs/spetrogram_hplot.png")

The power spectrum is easily computed and easily plotted. It is easy to set Welch's method, Barlett's method, 
or direct (no segmentation) PSD estimation. In this case we use Welch's method with 3 second windows and $0.5$ overlap.

```julia
psd = PSD(signal, eeg.fs, eeg.fs * 3,0.5)
plot(psd.freq, pow2db.(psd.spectrum), xlab="Frequency (Hz)", ylab="PSD (dB)", legend=false)
```

![Image]("imgs/psd.png")

To do: 

- [] Filtering
- [] Resampling

