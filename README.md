<img src="https://github.com/federicopirola/DBNcare/blob/main/notebooks/images/DBNcare.svg" alt="DBNcare package logo" width="227" height="227"/>
<h1>DBNcare: An R Package for Dynamic Bayesian Networks in Healthcare</h1>
<h2>Table of Contents</h2>
<ul>
  <li><a href="#introduction">Introduction</a></li>
  <li><a href="#installation">Installation</a></li>
  <li><a href="#usage">Usage</a></li>
  <li><a href="#contributing">Contributing</a></li>
  <li><a href="#license">License</a></li>
  <li><a href="#contact">Contact</a></li>
</ul>

## Introduction
<a name="introduction"></a>
<!-- Content for the Introduction section -->
Welcome to DBNcare, an evolving suite designed for Dynamic Bayesian Networks (DBNs). Recognizing the absence of a comprehensive toolset for DBNs, DBNcare is being developed to fill this gap. Our goal is to provide a full-featured library that supports the modeling, learning, inference, and visualization of dynamic processes. As DBNcare is still under active development, we continuously add new features and improvements. Whether conducting research or developing applications, DBNcare aims to become your ultimate solution for all things related to Dynamic Bayesian Networks.
## Installation
<a name="installation"></a>
<!-- Content for the Installation section -->
To install the DBNcare package, you will need to use R and the `devtools` package to install it directly from GitHub. Follow the steps below to get started:

1. Ensure you have R installed on your system. You can download it from [CRAN](https://cran.r-project.org/).
2. Install the `devtools` package if you haven't already:

```R
install.packages("devtools")
```

3. Use `devtools` to install DBNcare from GitHub:

```R
devtools::install_github("https://github.com/federicopirola/DBNcare")
```

After completing these steps, you'll be ready to use DBNcare for your Dynamic Bayesian Networks projects. As the package is still under development, make sure to check for updates regularly.
## Usage
<a name="usage"></a>
<!-- Content for the Usage section -->
The table below highlights the main features currently implemented in DBNcare. For more detailed information on these and additional features, please refer to the [documentation]() or explore this [notebook](https://github.com/federicopirola/DBNcare/blob/main/notebooks/dbn_main_features.Rmd).

| Feature             | Discrete | Continuous | Mixed     |
|---------------------|----------|------------|-----------|
| Missing Values      | ❌       | ❌         | ❌        |
| Parameter Learning  | ✔️       | ❌         | ❌        |
| Structure Learning  | ✔️       | ❌         | ❌        |
| Exact Inference     | ✔️       | ❌         | ❌        |
| Approximate Inference| ✔️      | ❌         | ❌        |
| Forecasting          | ✔️      | ❌         | ❌        |
| Sampling             | ✔️      | ❌         | ❌        |
## Contributing
<a name="contributing"></a>
<!-- Content for the Contributing section -->

## License
<a name="license"></a>
<!-- Content for the License section -->

## Contact
<a name="contact"></a>
<!-- Content for the Contact section -->
If you have any questions or feedback or need assistance with DBNcare, please feel free to reach out to us:

- 📧 Francesco Canonaco: [francesco.canonaco@minutia.ai](mailto:francesco.canonaco@minutia.ai)
- 📧 Federico Pirola: [federico.pirola@campus.unimib.it](mailto:federico.pirola@campus.unimib.it)
