## Overview

This project focusses on predicting whether nitrate concentration will be in the safe, warning or alert range, based on a variety of climate predictors, and the extent to which weather affects nitrate fluctuations. The results show a model within 60 and 80% accuracy based on the region, showing how models can often struggle to predict within smaller margins. This project including the geospatial analysis in it, can be used to both alert communities of high nitrate levels and advocate for future conservation efforts to reduce farm runoff. 

## Tools Used

- **RStudio/Python**: Data cleaning/transformation, ML model implementation, Geospatial analysis.
- **Excel**: Null value imputation.
- **Shiny**: App creation for interactivity.

## Data and Models Used

Data used includes the Ambient water testing data (available from the Iowa DNR) for the years 2022-2025 and daily weather data from the center for environmental information. Models used include polynomial regression and support vector regression, using monthly total rainfall and average temperature as predictors for the change in nitrate concentration response. The expected change was then added to the lag (previous month's concentration) and converted to the expected warning level. 

## Primary Findings

- **Inter-year Differences**: 2025 and 2024 had more tests above the EPA threshold of 10 mg/L $NO_3  - N$.
- **Seasonal Trends**: Nitrate concentration tended to peak in the late-spring/early-summer and leveled off by July/August.
- **Model Performance**: Neither the Polynomial nor SVR models performed better than the other.
- **Regional Performance**: Missouri river region had the highest F1 for the warning class. Larger waterways tended to have higher levels than smaller, regional waterways. 
