import imageio.v2 as imageio

#USER INPUTS_________________________________________________________________

path = r"C:\Users\nakaj\OneDrive\Desktop\2025Work\MassiveImpactTest\MergVisuals\Snapshots3\\"    
outputpath=r"C:\Users\nakaj\OneDrive\Desktop\2025Work\MassiveImpactTest\MergVisuals\Snapshots3\\"
prefix='3Dframe'
outputnumber1=0
outputnumber2=740
fps=17

#____________________________________________________________________________

def main():    # Main function generates the animation

    print('\nGenerating animation...\n')

    output = f'{outputpath}/Animation.mp4'
    writer = imageio.get_writer(output, fps=fps)

    for i in range(outputnumber1,outputnumber2+1):

        print(f'    Reading from {prefix}_{i:05d}.png')

        frame = imageio.imread(f'{path}/{prefix}_{i:05d}.png')
        writer.append_data(frame)

    writer.close()

    print(f'\n    Outputted Animation.mp4\n')
    


main()


