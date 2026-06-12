from pathlib import Path

path = Path('Data preparation classificaion TAG1-Copy1.ipynb')
text = path.read_text(encoding='utf-8')
replacements = {
    "dir_name='C:/Users/utraf/Desktop/16022022_RUSH MDCK EGFP_GPIFRa_50k_D5P5_500uMConsBiot/NATHALY/'": "dir_name='path/to/image_data/'",
    "filename = 'RUSH MDCK EGFP_GPIFRa_50k_D5P5_500uMConsBiot_60x1_t1.tif'": "filename = 'sample_image.tif'",
    "C:\\Users\\utraf\\AppData\\Local\\Temp\\ipykernel_20180\\3922955033.py:3: DeprecationWarning: <tifffile.imsave> is deprecated. Use tifffile.imwrite\\n": "C:\\Users\\<anon>\\AppData\\Local\\Temp\\kernel.py:3: DeprecationWarning: <tifffile.imsave> is deprecated. Use tifffile.imwrite\\n"
}
for old, new in replacements.items():
    if old not in text:
        print('Missing:', old)
path.write_text(text, encoding='utf-8')
print('Done.')
