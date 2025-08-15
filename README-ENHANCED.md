# pyKO Enhanced: Hybrid Spall & Interface Analysis

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**Enhanced version of pyKO with advanced spall and interface separation analysis capabilities**

## 🚀 What's New in This Version

### ✨ **Test 17: Hybrid Spall + Interface Separation Analysis**
- **Modular Analysis Framework:** Toggle spall, interface, FSV, and stress analysis independently
- **Dynamic Configuration:** Automatic YAML selection based on physics requirements
- **Advanced Visualization:** Custom pressure colormaps with physical meaning
- **Dual Spall Detection:** Both density-based and pressure-based spall identification
- **Comprehensive Parameter Display:** All simulation parameters extracted and displayed from YAML

### 🎯 **Key Features**
- **Custom Pressure Colormap:** Red (tension) → Gray (zero) → Blue (compression)
- **Enhanced Error Handling:** Robust data processing with safety checks
- **Detailed Configuration Guidelines:** Comprehensive YAML header documentation
- **Free Surface Velocity (FSV) Tracking:** Real-time impact surface monitoring
- **Stress Analysis:** Maximum compressive/tensile stress tracking
- **Interface Position Tracking:** Material boundary evolution

## 📊 **Analysis Capabilities**

### **Modular Analysis Toggles**
```python
# User-configurable analysis modules
ENABLE_SPALL_ANALYSIS = True        # Density + pressure-based spall detection
ENABLE_INTERFACE_ANALYSIS = True    # Material interface tracking  
ENABLE_FSV_ANALYSIS = True          # Free surface velocity analysis
ENABLE_STRESS_ANALYSIS = True       # Max stress analysis in target
```

### **Smart Configuration Selection**
- **With Interface Separation:** `test17-with-interface-separation.yml`
- **Without Interface Separation:** `test17-without-interface-separation.yml`
- **Automatic Selection:** Based on `ENABLE_INTERFACE_ANALYSIS` toggle

### **Advanced Visualization**
- **Eulerian x-t Diagrams:** Current position-based analysis
- **Lagrangian x-t Diagrams:** Initial position-based tracking
- **Custom Pressure Colormap:** Physically meaningful color representation
- **High-Resolution Plots:** Configurable DPI and time resolution

## 🔧 **Configuration Guidelines**

### **Mesh Selection**
```yaml
# Resolution targets:
# Shock physics:     2-5 μm/cell minimum
# Spall detection:   1-3 μm/cell (high resolution needed)
# Interface track:   2-4 μm/cell
# Computational:     5-10 μm/cell (faster, lower accuracy)

mat1:  # Example: 100 μm Al flyer
    mesh:
        cells  : 100          # 1 μm/cell (high resolution)
        xstart : -0.0001      # Start at -100 μm 
        length :  0.0001      # 100 μm thickness
```

### **Spall Parameters**
```yaml
# Material-specific fracture thresholds:
# Al:    ~200-500 MPa (2E8 - 5E8 Pa)
# Cu:    ~800-1200 MPa (8E8 - 1.2E9 Pa)
# Steel: ~1000-2000 MPa (1E9 - 2E9 Pa)

fracture:
    pfrac : 2.76E8        # Al spall threshold (276 MPa)
    nrhomin : 0.9         # Allow 10% density reduction
```

## 📁 **File Structure**

```
pyko-enhanced/
├── example-inputs/
│   ├── pko-test17-spall-interface.py           # Main hybrid analysis script
│   └── test17-spall-interface/
│       ├── test17-with-interface-separation.yml     # Full physics config
│       ├── test17-without-interface-separation.yml  # Spall-only config
│       └── test17-spall-interface.yml              # Base config
├── eos/                                        # Equation of state tables
├── pyko.py                                     # Core simulation engine
└── README-ENHANCED.md                          # This documentation
```

## 🏃‍♂️ **Quick Start**

### **1. Basic Impact Simulation**
```python
# Run Al flyer -> Cu target impact with full analysis
python pko-test17-spall-interface.py
```

### **2. Configure Analysis Modules**
```python
# Edit script header to enable/disable analysis modules
ENABLE_SPALL_ANALYSIS = True        # Spall detection
ENABLE_INTERFACE_ANALYSIS = True    # Interface tracking
ENABLE_FSV_ANALYSIS = True          # Free surface velocity
ENABLE_STRESS_ANALYSIS = True       # Stress analysis
```

### **3. Customize Material Properties**
```yaml
# Edit YAML files for different materials/configurations
# All parameters documented in file headers
```

## 📈 **Analysis Output**

### **Simulation Summary**
- Material configuration and properties
- Timing and resolution parameters
- Spall detection status (dual method)
- Maximum stress values and locations
- Interface separation detection

### **Visualization**
- Pressure x-t diagrams (Eulerian & Lagrangian)
- Free surface velocity vs. time
- Maximum stress evolution
- Material interface positions
- Density ratio evolution

## 🔬 **Physics Capabilities**

- **Elastic-Plastic Deformation:** Von Mises and hydrodynamic strength models
- **Shock Wave Propagation:** High-fidelity Lagrangian scheme
- **Spall Fracture:** Density and pressure-based detection
- **Interface Separation:** Material boundary tracking
- **Equation of State:** Mie-Gruneisen with tabular extensions

## 📚 **Original pyKO Information**

**Based on:** pyKO v0.6.1 by Sarah T. Stewart  
**Original Repository:** [https://github.com/ImpactsWiki/pyko](https://github.com/ImpactsWiki/pyko)  
**Documentation:** [https://impactswiki.github.io/pyko/](https://impactswiki.github.io/pyko/)  
**Citation:** Stewart, S. T. pyKO code v0.6.1, doi:10.5281/zenodo.8092348, 2023.

## 🤝 **Contributing**

Contributions welcome! This enhanced version builds upon the excellent foundation of the original pyKO code with additional analysis capabilities for shock physics research.

## 📄 **License**

GNU General Public License v3.0 - see [LICENSE](LICENSE) file for details.

## 🎓 **Citation**

If you use this enhanced version in your research, please cite both:
1. The original pyKO code: Stewart, S. T. pyKO code v0.6.1, doi:10.5281/zenodo.8092348, 2023.
2. This enhanced version: [Your repository citation when published]

---

**Enhanced by:** Piyush Wanchoo  
**Institution:** [Your Institution]  
**Contact:** [Your Email]
