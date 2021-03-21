# Display_PMA_Colormaps

Assume we already calculated the perfusion maps through some methods. Let’s call the generated perfusion maps as CBF and CBV (2D matrix).

**Requirements**: (everything is in the folder _/Display_PMA_Colormaps_) 
- Full PMA Color Lookup Table (CLT) from PMA software, _PMA_lut.csv_
- Function _select_colormap.m_  % select the PMA colormap that you want to display 
- Function _ctshow_pma.m_       % display images with defined colormap
- _CBF.mat, CBV.mat, mask.mat_  % required data for displaying (pre-calculated)

**Example**:
Run _example_ctshow_pma_colormap_full.m_

**Resulted image**:
![Example Output](https://github.com/yxiao009/Display_PMA_Colormaps/blob/master/results/result.jpg?raw=true)

**PMA**
[Source Page of PMA](http://asist.umin.jp/index-e.htm).
