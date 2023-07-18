# applied-statistics-project
Results: https://www.canva.com/design/DAFocS33n0o/hvD8xCp9VLIBE3RMownTUQ/edit?utm_content=DAFocS33n0o&utm_campaign=designshare&utm_medium=link2&utm_source=sharebutton

Ground Motion Models (GMMs) are complex equations used by seismologists to predict earthquake occurrence and the intensity of ground shaking. They play a crucial role in seismic hazard analysis and civil protection planning. In our university project, we focused on ITA18 GMM, the state-of-the-art model for predicting spectral acceleration in Italy.

Our goal was to develop a predictive model for estimating the Spectral Acceleration (SA) at unmeasured locations, considering spatial correlation and period dependence. Through geospatial plots, we observed a correlation between the magnitude of observed events and the style of faulting at different locations in Italy.

To analyze the influence of the style of faulting (SoF) on magnitude, we conducted ANOVA testing. The results indicated a significant relationship, leading us to implement three different models based on SoF data subdivision.

We also explored the relationship between SA and period (T) using contrast matrices. The findings revealed the importance of period as a crucial factor in our model.

Starting with a simple model, we gradually incorporated variables based on their significance, including SoF as a categorical regressor. The performance evaluation based on RMSE guided our selection of the optimal model.

We addressed collinearity issues by employing Ridge Regression with penalization adjustments for each period. However, the cross-validation consistently favored the reference model due to its physical robustness and interpretability.

By extending the ground motion models to a functional framework, we developed a robust and significant tool for predicting ground motion intensity across a continuous range of periods. This functional approach provides continuous estimation, inherent regularization, and implicit time-dependence.

Utilizing linear mixed effects models, we captured the internal variability among seismic events, which can be sporadic and distinct from one another.

In conclusion, our project demonstrated the enhancement of ground motion modeling and prediction by extending state-of-the-art models to a functional framework. We showcased the benefits of continuous estimation, regularization, and time-dependence. The inclusion of linear mixed effects models improved our understanding of internal variability among seismic events. This work contributes valuable insights for seismic analysis and prediction, supporting seismic hazard analysis and civil protection planning efforts in Italy.

Code and analysis by:
    erica.casassa@mail.polimi.it
    martina.fervari@mail.polimi.it
    mariagiovanna.latorraca@mail.polimi.it
    eliasraphael.roux@mail.polimi.it
    edoardo1.vitale@mail.polimi.it
