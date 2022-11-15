# Starpasta
A python implementation of Hurley et al. 2000, “Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity”
with the following modifications:
- Improved neutron star/black hole masses and electron-capture supernovae from Fryer et al. 2012, "COMPACT REMNANT MASS FUNCTION: DEPENDENCE ON THE EXPLOSION MECHANISM AND METALLICITY"
- Additional mass loss mechanisms from Belczynski et al. 2010, "ON THE MAXIMUM MASS OF STELLAR BLACK HOLES"
- Pair-instability supernovae from Belczynski et al. 2016, "The effect of pair-instability mass loss on black-hole mergers"
- Habitable Zone Boundaries from Kopparapu et al. 2014, "HABITABLE ZONES AROUND MAIN-SEQUENCE STARS: DEPENDENCE ON PLANETARY MASS"

## Use

Either starpasta.py or starpasta.exe; the former is a modifiable python script requiring a python installation with the numpy package, the latter is a standalone file created with the pyinstaller package.

Running either will open a command prompt window that will prompt you for the star's mass and metallicity.

Mass is in multiples of the sun's mass. The results should be reasonably accurate for the full evolution of stars between 0.8 and 150 solar masses. Main sequence evolution should be accurate for stars down to 0.5 solar masses and at least the early main sequence down to the minimum star mass of 0.08 solar masses, but post-MS evolution seems to deviate from expectations; you'll be prompted on whether you want to simulate full evolution or just the main sequence for these stars. These formulae have been used occasionally in some sources up to 300 solar masses, but accuracy isn't guaranteed and crashes may occur at low metallicity.

Metallicity is the portion of the star's mass composed of elements heavier than helium, which is roughly 0.02 for the sun; as with the original paper, the intended range is 0.0001 to 0.03.

The code should then run, print out a short summary of the results, and print out the full results to a .csv file; this shouldn't usually take more than a few seconds. Press 'enter' to close the window. 

The results .csv file, which will be named after the mass and metallicity, can be read using starpasta_out.xlsx: copy the contents of the .csv to the "Raw Data" tab of starpasta_out, and then you can use the "Results" section to examine it:

- The "Output Summary" box gives you a quick summary of the star’s evolution between stages and the duration of each one in millions of years.
- "Single Output" lets you specify a specific age of the star and see the parameters at that particular point in time. These values are all linearly interpolated between the 2 nearest timesteps, which should be a decent estimate because the Starpasta code tries to ensure that timesteps are always short enough that none of these values change too much (excepting special events like supernovae).
- The other boxes allow you to control graphs of the star’s whole evolution—both a Hertzprung–Russel diagram on the right and a user-defined graph below it, both of which display different stages in a star’s lifetime as different colors. Output Controls allows you to control which stages of the star’s lifetime are represented, or alternatively pick a start and end time to the data shown (this is necessarily rounded to the next timestep for the start time or the previous timestep for the end time). Starpasta will track the evolution of remnants (white dwarfs, neutron stars, black holes) for 10 trillion years, so you may need to exclude them or pick an appropriate end times to get reasonable graphs over time of the star’s evolution during its fusing lifetime.
- "Graph Control" below that will allow you to pick two parameters from the list, which will be charted against each other in the graph to its right. The graph will tend to automatically start its axes at 0, so if you want a more detailed graph you may have to adjust the settings.

## Notes

There are a number of ambiguities or oddities in the source code that I've done my best to reasonably interpret in my implementation:

- The main paper (Hurley et al. 2000) specifies that core mass during the HG should be checked against that in the previous timestep to ensure it isn’t unphysically shrinking due to mass lass, but then later says that unphysical core shrinkage should be avoided by holding the initial mass constant when necessary, which would seem to be redundant. I implement the latter solution, as it is discussed in a bit more length.
-	The core mass and stage timescale are interdependent during CHeB, so I take the last timestep’s core mass as an initial estimate for the current core mass, then iterate through calculating both values until there’s <0.1% variance between iterations.
-	Frankly there’s a lot of oddities in the CHeB stage that required various extra checks to prevent divide-by-zero errors in particular edge cases or account for the brief use of complex numbers, but I don't think any of these should change the code's normal behavior
-	The paper suggests that timestep lengths during TPAGB should be based on the difference between the current time and tinf2, but it seems in some cases that the TPAGB continues up to and past tinf2, with the result that the simulation gets trapped in an infinite loop at tinf2. I’m not sure if this is supposed to be possible but at any rate I set that particular timestep calculation to have a minimum of 1 year and the results in these cases seem to be reasonable.
-	For the small-envelope estimation, the paper says to estimate core mass during the GB as a zero-age naked helium star for M < MHeF and a white dwarf otherwise, but this is almost certainly a typo and should be the reverse.
-	Some low-mass, high metallicity stars can lose their envelopes on the EAGB and thus evolve to naked helium giants, but they immediately have core masses above the expected maximum mass for helium fusion to cease (which then throws errors during small envelope adjustment). I assume in these cases that they skip straight from EAGB to CO white dwarfs, even though the paper doesn’t mention this possibility.
-	Fryer et al. 2012 gives updates the post-supernova masses and accounting for electron-capture Sne, but is not explicit in how this would be implemented in the formulae, so I had to make some guesses: for stars with McBAGB between 1.83 and 2.25, I assume that the CO core converts completely to ONe and an EC Sne occurs if McCO reaches 1.38. All other stars supernova or collapse if they reach McSN, computed as in Hurley et al. 2000; those with lower McBAGB leave no remnant (I haven’t seen a parameter range for which this actually happens), those with higher are given a mass and type as described in Fryer et al. 2012, with McSN assumed to be the final CO core mass.
-	Kopparapu et al. 2014 gives HZ fits for stars with effective temperatures between 2600 and 7200 K; from their results, it’s clear we should expect the trends to higher Seff to continue for even hotter stars, but applying the given formula to much higher temperatures gives clearly unphysical results. Instead, I simply take the Seff values at 2600 and 7200 K and apply these for all lower and higher temperatures, respectively.
