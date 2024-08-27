# ACAP_JARA
Assessment of population trend for 9 priority populations of the Agreement for the Conservation of Albatrosses and Petrels (ACAP), analysed with the R package JARA (Just Another Red List Assessment) to estimate trends and likely Red List status inferred from population trends and population size.

This analysis was requested by Richard Phillips in August 2024 to support presentations at the seabird conference in Coimbra (September 2024) and for an ACAP meeting later in the year.

The analysis is based on count data and the the trend model consists of a state-space model for which both an observation process error and a state process error are separately defined. Priors for the observation process error were informed by the stated uncertainty associated with each count, and were set as follows:

Original:("High","H","high")= 0.25, ("Medio","Medium") = 0.4, ("unknown","Low", and any missing value) = 0.6
Modified:("High","H","high")= 0.1, ("Medio","Medium") = 0.2, ("unknown","Low", and any missing value) = 0.3

Any counts from sites with <2 counts were discarded because it is not possible to estimate a trend with a single point estimate. If a range was given for a certain count date and location (e.g. 200-400) then the the mean count of that range (e.g. 300) was used for analysis. The generation length was set to 20 years for procedural reasons, because ACAP requires a trend estimate over the past 20 years - care therefore needs to be taken when interpreting the IUCN Red List status assessments because for some species a generation length may be longer (e.g. Diomedea?) or shorter (e.g. Puffinus?) than the 20 years specified.
