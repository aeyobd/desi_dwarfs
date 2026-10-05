
files=(
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp00.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp01.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp02.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp03.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp04.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp05.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp06.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp07.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp08.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp09.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp10.fits'
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/catalogs/fastspec-iron-main-bright-nside1-hp11.fits'
)


for file in "${files[@]}"; do
  curl -O $file
done


