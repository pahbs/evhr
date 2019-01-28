#!/usr/bin/env python
#
# Get a list of attributes of a DSM from the corresponding XML file and calculate stereo angles
#
# Import and function definitions
import os, sys, math, osgeo, csv
from osgeo import ogr, osr, gdal

import datetime, time
from datetime import datetime
#gdal.AllRegister() #register all raster format drivers

# Function for calculating a 3x3 determinant
def det3(a1, b1, c1, a2, b2, c2, a3, b3, c3):
    res = a1*b2*c3+a2*b3*c1+a3*b1*c2-a1*b3*c2-a2*b1*c3-a3*b2*c1
    return res

def calc_stereoAngles(alpha1,theta1,alpha2,theta2,x1,y1,z1,x2,y2,z2,lat,lon):
    """
    alpha1  =   meanSatEl of image 1
    theta1  =   meanSatAz of image 1
    alpha2  =   meanSatEl image 2
    theta2  =   meanSatAz image 2
    x,y,z   =   satellite empheris

    References:
    -----------
    Jeong, J., C. Yang, and T. Kim. 2015. ?Geo-Positioning Accuracy Using Multiple-Satellite Images: IKONOS, Quickbird, and KOMPSAT-2 Stereo Images.?
    Remote Sensing 7 (4): 4549?4564.doi:10.3390/rs70404549.
    http://www.mdpi.com/2072-4292/7/4/4549

    Jeong, J., and T. Kim. 2015. ?Comparison of Positioning Accuracy of a Rigorous Sensor Model and Two Rational Function Models for Weak Stereo Geometry.?
    ISPRS Journal of Photogrammetry and Remote Sensing 108: 172?182. doi:10.1016/j.isprsjprs.2015.07.006.
    http://www.sciencedirect.com/science/article/pii/S0924271615001859

    Jeong, Jaehoon; Kim, Taejung. 2014. "Analysis of Dual-Sensor Stereo Geometry and Its Positioning Accuracy"
    Photogrammetric Engineering & Remote Sensing, Number 7 / July 2014, pp. 653-661(9)
    https://doi.org/10.14358/PERS.80.7.653

    Jeong, Jaehoon; Kim, Taejung. 2016. Quantitative Estimation and Validation of the Effects of the Convergence, Bisector Elevation, and Asymmetry Angles on the Positioning Accuracies of Satellite Stereo Pairs
    Photogrammetric Engineering & Remote Sensing Volume 82, Issue 8, August 2016, Pages 625-633
    https://www.sciencedirect.com/science/article/pii/S0099111216301021

    Jeong, Jaehoon. 2017. Use of weighted least squares for geo-positioning using dual-satellite image pairs
    International Journal of Remote Sensing Volume 38, 2017 - Issue 6
    http://www.tandfonline.com/doi/full/10.1080/01431161.2017.1285502

    Zhu, L., Umakawa, H., Guan, F., Tachibana, K., and Shimamura, H., 2008. Accuracy Investigation of Orthoimages Obtained from High Resolution Satellite Stereo Pairs,
    The International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, Beijing, Vol. XXXVII. PartB1. pp.1145-1148
    www.isprs.org/proceedings/XXXVII/congress/1_pdf/195.pdf

    Other sources and background:
    -----------
    http://www.geoimage.com.au/media/brochure_pdfs/DEMBrochure_FEB2015.pdf

    In the case of the VHR satellites with their pointable telescopes, the B/H ratio is not appropriate as a measure
    of the effectiveness of the stereo pair for DEM generation. In such cases, three angular measures of
    convergent stereo imaging geometry: the convergence angle, the asymmetry angle, and the bisector elevation angle (BIE) are used.

    These measure the geometrical relationship between two rays that intersect at a common ground point, one
    from the fore image and one from the aft image as shown in the diagram.

    Definitions:
    -----------
    Convergence Angle:
    The angle between two rays of a stereo pair
    The most important of the three stereo angles is the convergence and is the angle between the two rays in the
    convergence or epipolar plane. An angle between 30 and 60 degrees is ideal (<--- ideal for measuring what?? which heights? short trees??)

    Asymmetry Angle:
    The asymmetry angle is the angle between the bisector of the convergence angle and the projection of the local vertical onto the epipolar plane.
    Asymmetry describes the apparent offset from the centre view that a stereo pair has. For instance, a stereo pair
    with an asymmetry of 0? will have parallax due to elevations that appear equivalent in the left and right images.
    An asymmetrical collection is preferred as it gives a different look angle to discern ground features more
    accurately but should be under 20 deg.

    Bisector Elevation Angle:
    The obliqueness of the epipolar plane. BIE = 90 is orthogonal to ground surface
    The elevation angle of the bisector of the convergence angle
    The BIE angle is the angle between the horizontal plane and the epipolar plane and defines the amount of parallax that will
    appear in the vertical direction after alignment. The angle should be between 60 and 90 degrees.
    """

    # New Calculation
    # Using equations in J. Jeong and T Kim, 2016, Photogrammetric Engin. and RS, Vol. 82, No. 8

    # Degrees to radians
    dtr = math.atan(1.0) / 45.0

    # Set Earth Radius
    r = 6378137.   # WGS84 equatorial earth radius in meters

    # ECEF coordinate calculation
    # Position on ground, calc'd according to lat, lon, and earth radius
    # position is the center of the image
    x0 = r * math.cos(lat*dtr) * math.cos(lon*dtr)
    y0 = r * math.cos(lat*dtr) * math.sin(lon*dtr)
    z0 = r * math.sin(lat*dtr)

    # check the altitudes of the two satellite positions
    d1 = math.sqrt(x1*x1+y1*y1+z1*z1)
    d2 = math.sqrt(x2*x2+y2*y2+z2*z2)
    a = d1/d2
    x2 = x2 * a
    y2 = y2 * a
    z2 = z2 * a
    ##print '(x2,y2,z2) was adjusted by multiplying %s' %a

    # calculate theta1 and theta2
    # theta 1 is formed by two vectors U and V
    # U is from S1 to S2 and V is from S1 to the image center
    # S1 and S2 are positions of two cameras

    u1 = x1-x2
    u2 = y1-y2
    u3 = z1-z2
    du = math.sqrt(u1*u1 + u2*u2 + u3*u3)
    v1 = x1 - x0
    v2 = y1 - y0
    v3 = z1 - z0
    dv = math.sqrt(v1*v1 + v2*v2 + v3*v3)
    dd = (u1*v1 + u2*v2 + u3*v3) / (du * dv)
    theta1 = math.acos(dd) / dtr

    # theta2 is also formed by two vectors
    # U is from S2 to S1 and V is from S2 to the image center
    u1 = -u1
    u2 = -u2
    u3 = -u3

    v1 = x2 - x0
    v2 = y2 - y0
    v3 = z2 - z0
    dv = math.sqrt(v1*v1 + v2*v2 + v3*v3)
    dd = (u1*v1 + u2*v2 + u3*v3) / (du * dv)
    theta2 = math.acos(dd) / dtr
    ##print 'Theta1 and theta 2: %s and %s' %(theta1, theta2)

    # Convergence Angle (con_ang)
    con_ang = 180.0 - theta1 - theta2                                          # Equation 1 in the paper

    aa = math.sin(alpha1*dtr) * math.sin(alpha2*dtr)+ math.cos(alpha1*dtr) * math.cos(alpha2*dtr)* math.cos((theta1-theta2)*dtr)
    sc = (math.sin(alpha1*dtr) + math.sin(alpha2*dtr)) /(math.sqrt(2.0) * math.sqrt(1.0 + aa))

    # Bisector Elevation Angle (bie_ang)
    bie_ang = math.asin(sc)/dtr                                                # Equation 2 in the paper

    # Asymmetry Angle (asym_ang)
    asym_ang = (theta1 - theta2) * 0.5                                         # Equation 3 in the paper

    con_ang = round(con_ang, 3)
    bie_ang = round(bie_ang, 3)
    asym_ang = round(asym_ang, 3)

    ##printf, ou, format = '(6F10.4)', con_ang, asym_ang, bie_ang, con_ang2, asym_ang2, bie_ang2
    return (con_ang,bie_ang,asym_ang)


def main(imageDir):

    """
    imageDir        =   Top dir from which a pair of XMLs are read in

    Make lists of each set of XMLs for each catID in the iamgeDir name
    Calc angles for all combinations of XMLs from each catID list

    Output a CSV file with stereo angles and input XML info
    """
    # Create header
    hdr =   "IMG_L,IMG_R,"+\
            "YEAR,MONTH,DOY,"+\
            "TIME_L,TIME_R,"+\
            "SATEL_L,SATAZ_L,SUNEL_L,SUNAZ_L,ITVA_L,CTVA_L,ONVA_L,GSD_L,"+\
            "SATEL_R,SATAZ_R,SUNEL_R,SUNAZ_R,ITVA_R,CTVA_R,ONVA_R,GSD_R,"+\
            "CENTLON_L,CENTLAT_L,"+\
            "CENTLON_R,CENTLAT_R,"+\
            "EPHEMX_L,EPHEMY_L,EPHEMZ_L,"+\
            "EPHEMX_R,EPHEMY_R,EPHEMZ_R,"+\
            "ULLON_L,ULLAT_L,LLLON_L,LLLAT_L,URLON_L,URLAT_L,LRLON_L,LRLAT_L,"+\
            "ULLON_R,ULLAT_R,LLLON_R,LLLAT_R,URLON_R,URLAT_R,LRLON_R,LRLAT_R,"+\
            "ANG_CON,ANG_BIE,ANG_ASY\n"

    # Get pairname from input image dir
    baseDir, pairname = os.path.split(imageDir)
    print("\tDSM Info (epipolar geometry) for: %s" %(pairname))

    # Split pairname into catids
    cat1,cat2 = pairname.split("_")[2:] # get the last catIDs

    # Initialize lists
    cat1list = []
    cat2list = []

    # Create catid lists
    for root, dirs, files in os.walk(imageDir):
        for each in files:
            # Identify only xmls belonging to scenes
            if each.endswith('.xml') and 'P1BS' in each:
                if cat1 in each:
                    cat1list.append(each)
                if cat2 in each:
                    cat2list.append(each)

    # Name output csv with the pairname and put in output ASP dir
    outCSV = os.path.join(imageDir, "{}.csv".format(pairname))

    #--TMP NO CSV--#Open a CSV for writing
    #--TMP NO CSV--with open(outCSV,'wb') as csvfile:

        #--TMP NO CSV--# Write the header
        #--TMP NO CSV--#csvfile.write(hdr)
    #--TMP NO CSV--<INDENT EVERYTHING BELOW>
    # Get all combos of scenes from each catid strip:
    for leftXML in cat1list:
        for rightXML in cat2list:
            with open(os.path.realpath(os.path.join(imageDir,leftXML)), 'r') as file1, open(os.path.realpath(os.path.join(imageDir,rightXML)),'r') as file2:
                i = 0

                # Initialize vars - not necessarily stored for each stereopair (i.e., genYear, genMonth, genDOY will be the same for both)
                outline, catID = ('' for i in range(2))

                Year,Month,DOY,\
                tlc_dt,tlcTime,\
                meanSatEl,meanSatAz,meanSunEl,meanSunAz,meanITVA,meanCTVA,meanONVA,meanGSD,\
                ephemX,ephemY,ephemZ,\
                ullat,ullon,lllat,lllon,urlat,urlon,lrlat,lrlon,\
                maxLat,minLat,maxLon,minLon,centLat,centLon = (0 for i in range(30))

                catIDs, genDateCols, genTimeCols, SSGangles, ephemeris, centCoords, cornerCoords = ('' for i in range(7))

                # Loop through XML files
                for file in (file1,file2):
                    # Keep track of file
                    i += 1
                    # Read  XML line by line
                    for line in file.readlines():
                        ##print(line)
                        line2 = line.replace('<','>').split('>')[2]
                        # Get needed vars initialize above
                        if 'CATID' in line:
                            catID = str(line.replace('<','>').split('>')[2])
                        if 'TLCTIME' in line:
                            tlc_dt = datetime.strptime(str(line2).split('.')[0], '%Y-%m-%dT%H:%M:%S')

                            tlcTime = tlc_dt.strftime("%H:%M:%S")

                            Year = tlc_dt.year
                            Month = tlc_dt.month
                            DOY = tlc_dt.timetuple().tm_yday
                        if 'MEANSATEL' in line:
                            meanSatEl = float(line2)
                        if 'MEANSATAZ' in line:
                            meanSatAz = float(line2)
                        if 'MEANSUNEL' in line:
                            meanSunEl = float(line2)
                        if 'MEANSUNAZ' in line:
                            meanSunAz = float(line2)
                        if 'MEANINTRACKVIEWANGLE' in line:
                            meanITVA = float(line2)
                        if 'MEANCROSSTRACKVIEWANGLE' in line:
                            meanCTVA = float(line2)
                        if 'MEANOFFNADIRVIEWANGLE' in line:
                            meanONVA = float(line2)
                        if 'MEANPRODUCTGSD' in line:
                            meanGSD = float(line2)
                        else:
                            if 'MEANCOLLECTEDGSD' in line:
                                meanGSD = float(line2)
                        # Get Satellite Ephemeris using the first entry in EPHEMLISTList.
                        if '<EPHEMLIST>' in line and float(line.replace('<','>').replace('>', ' ').split(' ')[2]) == 1:
                            ephemX = float(line.replace('<','>').replace('>', ' ').split(' ')[3])
                            ephemY = float(line.replace('<','>').replace('>', ' ').split(' ')[4])
                            ephemZ = float(line.replace('<','>').replace('>', ' ').split(' ')[5])
                        if 'Upper Left' in line:
                            ullon = float(line.replace('(',')').split((')'))[1].split(',')[0])
                            ullat = float(line.replace('(',')').split((')'))[1].split(',')[1])
                        if 'Lower Left' in line:
                            lllon = float(line.replace('(',')').split((')'))[1].split(',')[0])
                            lllat = float(line.replace('(',')').split((')'))[1].split(',')[1])
                        if 'Upper Right' in line:
                            urlon = float(line.replace('(',')').split((')'))[1].split(',')[0])
                            urlat = float(line.replace('(',')').split((')'))[1].split(',')[1])
                        if 'Lower Right' in line:
                            lrlon = float(line.replace('(',')').split((')'))[1].split(',')[0])
                            lrlat = float(line.replace('(',')').split((')'))[1].split(',')[1])
                        if 'ULLON' in line:
                            ullon = float(line2)
                        if 'ULLAT' in line:
                            ullat = float(line2)
                        if 'URLON' in line:
                            urlon = float(line2)
                        if 'URLAT' in line:
                            urlat = float(line2)
                        if 'LLLON' in line:
                            lllon = float(line2)
                        if 'LLLAT' in line:
                            lllat = float(line2)
                        if 'LRLON' in line:
                            lrlon = float(line2)
                        if 'LRLAT' in line:
                            lrlat = float(line2)

                        maxLat = max(ullat,urlat,lllat,lrlat)
                        minLat = min(ullat,urlat,lllat,lrlat)
                        maxLon = max(ullon,urlon,lllon,lrlon)
                        minLon = min(ullon,urlon,lllon,lrlon)
                        centLat = minLat + (maxLat - minLat)/2
                        centLon = minLon + (maxLon - minLon)/2

                    # If only 1 catID, copy these vars for the first file.
                    if i == 1:
                        meanSatAz_1 = meanSatAz
                        meanSatEl_1 = meanSatEl
                        ephemX_1 = ephemX
                        ephemY_1 = ephemY
                        ephemZ_1 = ephemZ

                    # From both images:
                    # Get scene names instead of catids
                    Names        = leftXML + ',' + rightXML + ','

                    # gather date (not needed for both, so dont use += )
                    genDateCols  = str(Year) + ',' + str(Month) + ',' + str(DOY) + ','
                    genTimeCols += tlcTime + ','

                    # gather Sun-Sensor Geometry Angles
                    SSGangles   +=  str(meanSatEl)  + ',' + str(meanSatAz) + ',' + \
                                    str(meanSunEl)  + ',' + str(meanSunAz) + ',' + \
                                    str(meanITVA)  + ',' + str(meanCTVA)  + ',' + str(meanONVA) + ',' + str(meanGSD) + ','

                    # gather satellite ephemeris XYZ
                    ephemeris   +=  str(ephemX)     + ',' + str(ephemY) + ',' + str(ephemZ) + ','

                    # gather center coords
                    centCoords  +=  str(centLon)    + ',' + str(centLat) + ','

                    # gather corners
                    cornerCoords+=  str(ullon)      + ',' + str(ullat) + ',' + \
                                    str(lllon)      + ',' + str(lllat) + ',' + \
                                    str(urlon)      + ',' + str(urlat) + ',' + \
                                    str(lrlon)      + ',' + str(lrlat) + ','

                # Calc stereo angles
                stereoAngs = calc_stereoAngles(meanSatEl_1,meanSatAz_1,meanSatEl,meanSatAz,ephemX_1,ephemY_1,ephemZ_1,ephemX,ephemY,ephemZ,centLat,centLon)

                out_attribute_line = Names + genDateCols + genTimeCols + SSGangles + centCoords + ephemeris + cornerCoords + str(stereoAngs[0]) + "," + str(stereoAngs[1]) + "," + str(stereoAngs[2])+'\n'

                try:
                    with open(outCSV,'wb') as csvfile:
                        # Write the header
                        csvfile.write(hdr)
                        # Write attributes
                        csvfile.write(out_attribute_line)
                except Exception, e:
                    print "\tFailed to write CSV of dsm info: %s" %outCSV

                # Write line

                #--TMP NO CSV--csvfile.write(out_attribute_line)

                #--TMP NO CSV--print("\tOuput CSV file: %s" %(outCSV))
                print("\tConvergence Angle          = " + str(stereoAngs[0]))
                print("\tBisector Elevation Angle   = " + str(stereoAngs[1]))
                print("\tAsymmetry Angle            = " + str(stereoAngs[2]))

                return(str(stereoAngs[0]), str(stereoAngs[1]), str(stereoAngs[2]), hdr, out_attribute_line)


if __name__ == '__main__':
    main()