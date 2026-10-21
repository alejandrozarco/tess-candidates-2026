"""Extract SPOC LC FITS -> npz (own code). Usage: v1_extract.py <dl_dir> <out_dir>"""
import sys, os, glob, numpy as np
from astropy.io import fits
src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for f in sorted(glob.glob(os.path.join(src,'*lc.fits'))):
    fn=os.path.basename(f); fast='fast' in fn
    with fits.open(f, memmap=False) as h:
        hd0=h[0].header; hd=h[1].header; d=h[1].data
        sec=hd0['SECTOR']; tic=hd0['TICID']
        t=np.array(d['TIME'],float); q=np.array(d['QUALITY'],int)
        pf=np.array(d['PDCSAP_FLUX'],float); pe=np.array(d['PDCSAP_FLUX_ERR'],float)
        sf=np.array(d['SAP_FLUX'],float); se=np.array(d['SAP_FLUX_ERR'],float); bk=np.array(d['SAP_BKG'],float)
        meta=dict(sector=sec, tic=tic, fast=fast, crowd=hd.get('CROWDSAP',np.nan), flfr=hd.get('FLFRCSAP',np.nan),
                  procver=hd0.get('PROCVER',''), dataver=hd0.get('DATA_REL',-1), exptime=hd.get('TIMEDEL',np.nan)*86400, cam=hd0.get('CAMERA'), ccd=hd0.get('CCD'))
    good=np.isfinite(t)&np.isfinite(pf)&np.isfinite(sf)&(q==0)
    tag='%d_s%03d_%s'%(tic,sec,'20s' if fast else '2m')
    np.savez(os.path.join(out,tag+'.npz'), t=t[good], pdc=pf[good]/np.nanmedian(pf[good]), pdce=pe[good]/np.nanmedian(pf[good]),
             sap=sf[good]/np.nanmedian(sf[good]), sape=se[good]/np.nanmedian(sf[good]), bkg=bk[good]/np.nanmedian(sf[good]), meta=np.array([meta],dtype=object))
    print(tag, good.sum(), '%.3f %.3f'%(meta['crowd'],meta['flfr']), meta['procver'], meta['dataver'], round(meta['exptime'],1), 'rms_pdc_ppm=%.0f'%(1e6*np.nanstd(np.diff(pf[good]/np.nanmedian(pf[good])))/np.sqrt(2)))
