// Akima interpolation helper function
function calcA(i, xy) {
  if (xy.length < 4) return 0;
  
  const n = xy.length;
  if (i <= 1 || i >= n - 2) return 0;
  
  const m1 = (xy[i].efficiency - xy[i-1].efficiency) / (xy[i].energy - xy[i-1].energy);
  const m2 = (xy[i+1].efficiency - xy[i].efficiency) / (xy[i+1].energy - xy[i].energy);
  const m3 = (xy[i+2].efficiency - xy[i+1].efficiency) / (xy[i+2].energy - xy[i+1].energy);
  const m0 = i >= 2 ? (xy[i-1].efficiency - xy[i-2].efficiency) / (xy[i-1].energy - xy[i-2].energy) : m1;
  
  const dm1 = Math.abs(m2 - m1);
  const dm2 = Math.abs(m0 - m1);
  const dm3 = Math.abs(m3 - m2);
  
  if (dm1 + dm3 === 0) {
    return (m1 + m2) / 2;
  }
  
  return (dm3 * m1 + dm1 * m2) / (dm1 + dm3);
}

// Akima interpolation function
function akimaInterpolate(z, xy) {
  if (!xy || xy.length === 0) return 0;
  if (xy.length === 1) return xy[0].efficiency;
  
  // Find the interval
  let i = 0;
  for (i = 0; i < xy.length - 1; i++) {
    if (z <= xy[i + 1].energy) break;
  }
  
  if (i >= xy.length - 1) return xy[xy.length - 1].efficiency;
  if (i === 0 && z < xy[0].energy) return xy[0].efficiency;
  
  const x0 = xy[i].energy;
  const x1 = xy[i + 1].energy;
  const y0 = xy[i].efficiency;
  const y1 = xy[i + 1].efficiency;
  
  const h = x1 - x0;
  if (h === 0) return y0;
  
  const t = (z - x0) / h;
  
  const a0 = calcA(i, xy);
  const a1 = calcA(i + 1, xy);
  
  return y0 + t * (a0 * h + t * (3 * (y1 - y0) - (2 * a0 + a1) * h + t * ((a0 + a1) * h - 2 * (y1 - y0))));
}

// Exponential of log power series efficiency calculation
function expOfLogPowerSeriesEfficiency(energy, coefs) {
  if (!coefs || coefs.length === 0) return 0;
  if (energy <= 0) return 0;
  
  const logE = Math.log(energy);
  let sum = 0;
  let logEPower = 1;
  
  for (let i = 0; i < coefs.length; i++) {
    sum += coefs[i] * logEPower;
    logEPower *= logE;
  }
  
  return Math.exp(sum);
}

// FWHM calculation function
function peakResolutionFWHM(energy, form, pars) {
  if (!pars || pars.length === 0) return 0;
  
  switch (form) {
    case 'kGadrasResolutionFcn':
    {
      if (pars.length < 3) return 0;
      
      const A = pars[0];
      const B = pars[1]; 
      const C = pars.length > 2 ? pars[2] : 0;
      
      if (energy <= 0) return 0;
      
      const sqrtE = Math.sqrt(energy);
      const val = A + B * sqrtE + C * energy;
      
      if (val <= 0) return 0;
      return Math.sqrt(val);
    }

    case 'kSqrtEnergyPlusInverse':
    {
      if (pars.length < 3) return 0;
      if (energy <= 0) return 0;
      
      const val2 = pars[0] + pars[1] * energy + pars[2] / energy;
      return val2 > 0 ? Math.sqrt(val2) : 0;
    }

    case 'kConstantPlusSqrtEnergy':
    {
      if (pars.length < 2) return 0;
      if (energy <= 0) return 0;
      
      return pars[0] + pars[1] * Math.sqrt(energy);
    }

    case 'kSqrtPolynomial':
    {
      if (pars.length === 0) return 0;
      if (energy <= 0) return 0;
      
      // Use Horner's method to evaluate polynomial
      let val = pars[pars.length - 1];
      for (let i = pars.length - 2; i >= 0; i--) {
        val = val * energy + pars[i];
      }
      
      return Math.sqrt(val);
    }

    default:
      throw new Error('Unknown FWHM form: ' + form);
  }
}

// DetectorPeakResponseJS class for handling efficiency and FWHM calculations
class DetectorPeakResponseJS {
  constructor(data) {
    this.data = data;
  }
  
  // Check if efficiency data is available
  hasEfficiency() {
    return this.data && this.data.efficiency && this.data.efficiency.form;
  }
  
  // Check if FWHM data is available
  hasFwhm() {
    return this.data && this.data.fwhm && this.data.fwhm.form;
  }
  
  // Get energy extent
  getEnergyExtent() {
    if (this.data && this.data.extent) {
      return this.data.extent;
    }
    return [50, 3000]; // Default fallback
  }
  
  // Calculate efficiency at given energy
  efficiency(energy) {
    if (!this.hasEfficiency()) {
      return null;
    }
    
    const effData = this.data.efficiency;
    
    if (effData.form === 'kEnergyEfficiencyPairs') {
      const energyInUnits = energy / effData.energyUnits;
      return akimaInterpolate(energyInUnits, effData.pairs);
    } else if (effData.form === 'kExpOfLogPowerSeries') {
      const energyInUnits = energy / effData.energyUnits;
      return expOfLogPowerSeriesEfficiency(energyInUnits, effData.coefficients);
    }
    
    return null;
  }
  
  // Calculate FWHM at given energy
  fwhm(energy) {
    if (!this.hasFwhm()) {
      return null;
    }
    
    try {
      return peakResolutionFWHM(energy, this.data.fwhm.form, this.data.fwhm.coefficients);
    } catch (error) {
      return null;
    }
  }
}

// Make the class available globally
if (typeof window !== 'undefined') {
  window.DetectorPeakResponseJS = DetectorPeakResponseJS;
}

// For Node.js environments
if (typeof module !== 'undefined' && module.exports) {
  module.exports = DetectorPeakResponseJS;
}
