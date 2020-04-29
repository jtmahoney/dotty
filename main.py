import argparse
import xml.etree.ElementTree as et
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
import matplotlib.gridspec as gridspec
from skimage import io
import skimage
import numpy as np
import pandas as pd



def get_xml_dots(single_ch_image, xml_file_path, marker_names, box_size):
    tree = et.parse(xml_file_path)
    root = tree.getroot()

    column_names = ['x_coord', 'y_coord', 'x_min', 'x_max', 'y_min', 'y_max', 'box_size', 'marker_type', 'intensity_value']
    df_table = pd.DataFrame(columns=column_names)

    for mtype in root[1][1:]:
        marker_number = int(mtype.find('Type').text)
        for mark in mtype.findall('Marker'):
            if mark is not None:
                xs = int(mark.find('MarkerX').text)
                ys = int(mark.find('MarkerY').text)
                x_min = int(xs-(box_size/2))
                x_max = int(xs+(box_size/2))
                y_min = int(ys-(box_size/2))
                y_max = int(ys+(box_size/2))
                crop_box = single_ch_image[y_min:y_max, x_min:x_max]
                intensity_value = int(np.average(crop_box))
                df_table = df_table.append({'x_coord':xs,
                                            'y_coord':ys,
                                            'x_min':x_min,
                                            'x_max':x_max,
                                            'y_min':y_min,
                                            'y_max':y_max,
                                            'box_size':box_size,
                                            'marker_type':marker_names[marker_number-1],
                                            'intensity_value':intensity_value},
                                            ignore_index=True)
    return(df_table)

def plot_hist_img(data_table,color_dictionary,limits,full_image,display_channel,marker_size):
    df=data_table
    markers=list(df.marker_type.unique())
    df = df[(df['intensity_value'] >= limits[0]) & (df['intensity_value'] <= limits[1])]
    print('Number of data points: ', df.shape[0])
    gs = gridspec.GridSpec(3,3)
    fig = plt.figure(figsize=(30,40),tight_layout=True)

    ### PLOT THE HISTOGRAM ###
    ax0 = fig.add_subplot(gs[0,0])
    ax0.set_title('Intensity Histogram (filtered to n=%s)'%df.shape[0], fontsize=20)
    ax0.set_ylabel('n cells/bin', fontsize=10)
    ax0.set_xlabel('ave. pixel value', fontsize=10)
    ax0.hist(df['intensity_value'], bins=int(((limits[1]-limits[0])/20)))
    ax0.set_xlim(0,255)
    
    ### PLOT THE IMAGE ###
    ax1 = fig.add_subplot(gs[0:,1:])
    ax1.set_title('Dot Overlay Image (filtered to n=%s)'%df.shape[0], fontsize=20)
    ax1.imshow(full_image[display_channel])
    for index,row in df.iterrows():
        c = color_dictionary[row['marker_type']]
        x = row['x_coord']
        y = row['y_coord']
        circ = Circle((x,y),marker_size,color=c)
        ax1.add_patch(circ)
    
    ### PLOT THE GRAPH ###
    marker_array = np.zeros((len(markers),255), dtype=int)
    for mr in range(len(markers)):
        dfmark=data_table[data_table.marker_type == markers[mr]]
        for v in range(255):
            count_over_v=dfmark[dfmark.intensity_value >= v]
            marker_array[mr,v]=len(count_over_v.index)
    axCo = fig.add_subplot(gs[1,0])
    axCo.set_title('Cell Type Intensity (all)', fontsize=20)
    axCo.set_ylabel('n cells/subtype', fontsize=10)
    axCo.set_xlabel('ave. pixel value', fontsize=10)
    for mr in range(len(markers)):
        colorCo = color_dictionary[markers[mr]]
        axCo.plot(marker_array[mr,:], colorCo)
    
    ### PLOT THE KEY/LEGEND ###
    k = df['marker_type'].value_counts()
    axKey = fig.add_subplot(gs[2,0])
    axKey.set_title('ColorKey', fontsize=20)
    axKey.set_xlim(0,3)
    axKey.set_ylim(0,k.shape[0]+1)
    axKey.xaxis.set_visible(False)
    axKey.yaxis.set_visible(False)
    ypos=k.shape[0]
    for ke, va in k.items():
        kcirc = Circle((1,ypos),0.2,color=color_dictionary[ke])
        axKey.add_patch(kcirc)
        text_string = '{} ( {} )'.format(ke,va)
        ktext = axKey.text(2,ypos,text_string,fontsize=18)
        ypos=ypos-1
    
    return(fig)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    ### ANALYZE IMAGE USING XML COUNTS ###
    parser.add_argument('-x', '--xml_path', required=True, type=str, help='Path to CellCounter .xml file')
    parser.add_argument('-i', '--image_path', required=True, type=str, help='Path to multichannel .tif file')
    parser.add_argument('-m', '--marker_names', required=True, nargs='+', help='Names for each marker type, separated by a space')
    parser.add_argument('-s', '--syfp_channel', default=1, type=int, help='Channel number for vector (zero indexed)')
    parser.add_argument('-b', '--box_size', default=3, type=int, help='Size of grid used to measure vector expression level')
    parser.add_argument('-pkl', '--save_pickle', default=False, type=bool, help='True-False: save mesurements to .pkl file')

    ### MAKING A FIGURE ###
    parser.add_argument('-fig', '--make_fig', default=True, type=bool, help='Would you like to make a figure? True-False')
    parser.add_argument('-c', '--colors', default=['#2bfa05', '#fb04e9', '#ecf807', '#e32400', '#fec700'], nargs='+', help='Color for each marker separated by a space')
    parser.add_argument('-dc', '--display_channel', default=1, type=int, help='Channel# to display with markers on top (zero indexed)')
    parser.add_argument('-tl', '--threshold_low', default=0, type=int, help='Set the threshold for dimmest vector')
    parser.add_argument('-th', '--threshold_high', default=255, type=int, help='Set threshold for brightest vector')
    parser.add_argument('-dot', '--fig_dot_size', default=10, type=int, help='Size of dot to display on image')
    parser.add_argument('-save', '--save_fig_path', default='./figure.png', type=str, help='Path to save the figure')


    args = parser.parse_args()

    image = io.imread(args.image_path)
    fluor_img = image[args.syfp_channel]
    dtable = get_xml_dots(fluor_img, args.xml_path, args.marker_names, args.box_size)

    if args.save_pickle:
        dtable.to_pickle('datatable.pkl')

    if args.make_fig:
        color_dic = dict(zip(args.marker_names, args.colors))
        figgy = plot_hist_img(data_table=dtable,
                        color_dictionary=color_dic,
                        limits=(args.threshold_low, args.threshold_high),
                        full_image=image,
                        display_channel=args.display_channel,
                        marker_size=args.fig_dot_size)
        figgy.savefig(args.save_fig_path)
        
