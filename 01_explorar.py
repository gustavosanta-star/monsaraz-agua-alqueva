import ee
import urllib.request

ee.Initialize(project="projetobeta-513108")

aoi = ee.Geometry.Rectangle([-7.43, 38.40, -7.33, 38.49])

col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
       .filterBounds(aoi)
       .filterDate("2024-06-01", "2024-09-30")
       .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 10))
       .sort("CLOUDY_PIXEL_PERCENTAGE"))

print("Imagens encontradas:", col.size().getInfo())

img = col.first()
data = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd").getInfo()
print("Imagem escolhida, data:", data)

url = img.getThumbURL({
    "bands": ["B4", "B3", "B2"],
    "min": 0, "max": 3000,
    "region": aoi, "dimensions": 1024, "format": "png",
})

urllib.request.urlretrieve(url, "monsaraz_rgb.png")
print("Guardado: monsaraz_rgb.png")