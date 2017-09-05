from astropy.io import fits
from matplotlib import pyplot as p
import numpy as n
from scipy import integrate
from scipy import optimize
from scipy.special import wofz

def int_line_flux(x, y, x1, x2):
  """
  PURPOSE: Calculate the integrated area (using the trapezoidal rule) 
  		 of a spectral line while subtracting the continuum level

  ARGUMENTS:
    x, y: numpy ndarrays with (hopefully) matching sizes
    x1: lower limit of integral
    x2: upper limit of integral

  RETURNS: integrated line flux
  """
  #pick out the new x and y
  idx = (x1 <= x) & (x <= x2)
  x,y = x[idx], y[idx]
  #subtract out the continuum (the mean value between the bounds)
  y -= ((y[0]+y[(y.size-1)])/2.0)
  #could use scipy integrate at this point
  #print "simpsons rule for lambda = %i : %.3e" % (n.mean(x),integrate.simps(y,x))
  #or use the trapezoidal rule
  dx = n.diff(x)
  dy = n.diff(y)
  #take off the last element of y
  y = y[:-1]
  #use the trapezoidal rule
  area = n.sum(.5*dx*dy + dx*y)
  return area

def electron_temp(I4363,I4959,I5007):
  """
    PURPOSE: Calculate the electron temperature in the HII region
             using OIII line ratios

    ARGUMENTS:
      I4363,I4959,I5007: integrated flux of the OIII lines

    RETURNS: electron temperature in the HII region
  """
  flux_ratio = I4363/(I4959+I5007)
  print "OIII flux ratio = ", flux_ratio
  Te = -33000/(n.log(flux_ratio) - n.log(0.14))
  return Te

def gauss_fit(x,y,x1,x2):
  """
  PURPOSE: Fit a gaussian curve to data

  ARGUMENTS:
    x, y: numpy ndarrays with (hopefully) matching sizes
    x1: lower limit of curve fit
    x2: upper limit of curve fit

  RETURNS: best fit parameters
  """
  #pick out the new x and y
  idx = (x1 <= x) & (x <= x2)
  x,y = x[idx], y[idx]

  gaussfunc = lambda x,a,x0,sigma: a*n.exp(-(x-x0)**2/(2*sigma**2))
  #p0 = [a,x0,sigma]
  popt,pcov = optimize.curve_fit(gaussfunc,x,y,p0=[n.max(y),n.mean(x),(.5*(x.max()-x.min()))])
  
  print "optimized Gaussian fit for lambda = %i : " % n.mean(x), popt
  
  #overplot the fit with the data
  xfit = n.arange(x1,x2,0.1)
  yfit = n.copy(xfit)*0.0
  yfit = popt[0]*n.exp(-(xfit-popt[1])**2/(2*popt[2]**2))
  p.plot(x,y)
  p.plot(xfit,yfit)
  p.show()
  return popt

def lorentzian_fit(x,y,x1,x2):
  """
  PURPOSE: Fit a Lorentian profile to data

  ARGUMENTS:
    x, y: numpy ndarrays with (hopefully) matching sizes
    x1: lower limit of curve fit
    x2: upper limit of curve fit

  RETURNS: best fit parameters
  """
  #pick out the new x and y
  idx = (x1 <= x) & (x <= x2)
  x,y = x[idx], y[idx]
  y = y/y.mean()
  lorentzfunc = lambda x,x0,gamma: gamma/(2.*n.pi*((x-x0)**2 + (.5*gamma)**2))
  #popt = [x0,gamma]
  popt,pcov = optimize.curve_fit(lorentzfunc,x,y,p0=[n.mean(x),(.8*(x.max()-x.min()))])
  
  print "optimized Lorentian fit for lambda = %i : " % n.mean(x), popt
  
  #overplot the fit with the data
  xfit = n.arange(x1,x2,0.1)
  yfit = n.copy(xfit)*0.0
  yfit = popt[1]/(2.*n.pi*((xfit-popt[0])**2 + (.5*popt[1])**2))
  y = y*y.mean()
  p.plot(x,y)
  p.plot(xfit,yfit)
  p.show()
  return popt

def voigt(x, *p):
    A, mu, sigma, gamma = p
    z = ((x-mu) + 1j*gamma) / (sigma *n.sqrt(2.0))
    return A * n.real(wofz(z))

def voigt_fit(x,y,x1,x2):
  """
  PURPOSE: Fit a Voigt profile to data

  ARGUMENTS:
    x, y: numpy ndarrays with (hopefully) matching sizes
    x1: lower limit of curve fit
    x2: upper limit of curve fit

  RETURNS: best fit parameters
  """
  #pick out the new x and y
  idx = (x1 <= x) & (x <= x2)
  x,y = x[idx], y[idx]
  y /= y.mean()
  #voigtfunc = lambda x, A, x0, gamma, sigma: A * n.real(wofz(((x-x0) + 1j*gamma) / (sigma *n.sqrt(2.0))))
  #popt = [A,x0,gamma,sigma]
  popt,pcov = optimize.curve_fit(voigt,x,y,p0=[y.max(),x.mean(),1.0, 1.0])
  
  print "optimized Lorentian fit for lambda = %i : " % n.mean(x), popt
  
  #overplot the fit with the data
  xfit = n.arange(x1,x2,0.1)
  #yfit = n.copy(xfit)*0.0
  yfit = voigt(xfit,popt[0],popt[1],popt[2],popt[3])
  y *= y.mean()
  p.plot(x,y)
  p.plot(xfit,yfit)
  p.show()
  return popt

def main():
  #open the data file and get header parts
  spec_file = 'hIIspec.fits'
  hdulist = fits.open(spec_file)
  spec = hdulist[0].data
  crdelt = hdulist[0].header['cdelt1']
  crval = hdulist[0].header['crval1']
  hdulist.close()
  
  #make a matching wavelength array
  waves = n.arange(0,spec.size)*crdelt+crval
  #plot the spectrum
  p.plot(waves,spec)
  p.title('Some spectrum plot')
  p.xlabel('Wavelength (angstroms)')
  p.ylabel('Intensity')
  #save the plot
  p.savefig('spectrum')
  p.close()
  
  #integrate the OIII emission lines
  f4363 = int_line_flux(waves,spec,4361,4379)
  f4959 = int_line_flux(waves,spec,4957,4976)
  f5007 = int_line_flux(waves,spec,5005,5025)
  #calculate the electron temperature
  print "Te = %.4e" % electron_temp(f4363,f4959,f5007)

  #integrate the SII emission lines
  f6716 = int_line_flux(waves,spec,6717,6735)
  f6731 = int_line_flux(waves,spec,6734,6750)
  print "SII flux ratio: %.3e" % (f6716/f6731)
  #fit a Gaussian to the lines
  #gauss_fit(waves,spec,4361,4379)
  #gauss_fit(waves,spec,4957,4976)
  #gauss_fit(waves,spec,5005,5025)
  #gauss_fit(waves,spec,6717,6735)
  #gauss_fit(waves,spec,6734,6750)
  #fit a Lorentzian
  #lorentzian_fit(waves,spec,4361,4379)
  #lorentzian_fit(waves,spec,4957,4976)
  #lorentzian_fit(waves,spec,5005,5025)
  #lorentzian_fit(waves,spec,6717,6735)
  #lorentzian_fit(waves,spec,6734,6750)
  #fit a Voigt
  #voigt_fit(waves,spec,4361,4379)
  #voigt_fit(waves,spec,4957,4976)
  #voigt_fit(waves,spec,5005,5025)
  #voigt_fit(waves,spec,6717,6735)
  #voigt_fit(waves,spec,6734,6750)

if __name__ == '__main__':
    main()

