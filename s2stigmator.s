// $BACKGROUND$
//Astigmatism & Magnification distortion Correction

Number GMS = 0 //a variable to mark the GMS version, it is determined in main function. In GMS 1.x, size of images/array must be the power of 2.
Number mainImageID = 0 //a variable to recaord the image ID
Number voltage = 300
Number cs = 2
Number boxsize = 512  //boxsize for grid-boxing
Number clip = 512  //clip a smaller central region from the single ring image
Number oversample = 256 //In GMS 1.x, size of images/array/oversample must be the power of 2.
Number length = boxsize*oversample 
//Number oversample = length/boxsize
//Number oversample = 1
Number converter = 1 //DM3 saves the apix in nm, MRC saves the apix in um, convert it to A
Number apix = 1
Number highResolution = 5.5   //most microscopes may have the highest resolution 3~4A. Should not simply set it to Nyquist
Number lowX, highX
Number wavelen = 12.2639/sqrt(voltage * 1000.0 + 0.97845 * voltage * voltage)
String imageName
Number totalAstig = 100 //in nm
Number labelPosX = 0, labelPosY = 0, removePeakLabel = 0, removeComponents = 0
//Number fineFactor = 10
Number fineFactor2 = 1 //for dst20, if 1, the image of polar coordinates has 1*360 columes
Number fineFactor3 = 10 //for dst30, finer sampling than dst20
Number numWdgFactor = 1
Number numWdg1 = 256*numWdgFactor //The number of wedges when sampling on the single ring. It is equal to the number of points we sample on the single ring
Number verbose = 1  //0: no display, 1: display 2 final plots, 2: display all plots and figures

//Mark the ctffind3 results on the target plot
Number astigCtffind3 = 0   //in nm
Number angCtffind3 = 0    //in degree
Number s4Defocus = 0 //in nm
Number correctedDefocus = 0  //in nm

//Save result to a txt file
Number wMeanDf, wAstig, wAngle
Number fileID
Number saveResults = 0
Number nHistory = 500
Image historyX := realImage("historyX", 4, nHistory, 1)
Image historyY := realImage("historyY", 4, nHistory, 1)
Number bestX = 0, bestY = 0
Number size = 500, distance2Origin = size

//grayscale in the points
Number nGrayScale = 0
Number grayScaleLow = 0 //0 is black
Number grayScaleHigh = 1 


//other important parameters for recording
Number beamShiftX, beamShiftY //EMGetBeamShift
Number calBeamShiftX, calBeamShiftY //EMGetCalibratedBeamShift
Number beamTiltX, beamTiltY //EMGetBeamTilt
Number calBeamTiltX, calBeamTiltY //EMGetCalibratedBeamTilt
Number brightness //EMGetBrightness
Number getMag //EMGetMagnification
Number getFocus, getCalFocus //EMGetFocus, EMGetCalibratedFocus
Number condStigX, condStigY //EMGetCondensorStigmation
Number objStigX, objStigY //EMGetObjectiveStigmation
Number calObjStigX, calObjStigY //EMGetCalibratedObjectiveStigmation
Number spotSize //EMGetSpotSize





linePlotImageDisplay lpid

class idcListenFFT
  {
	  Number num_chg, img_id
	  //String task
	  Image imgCloned, imgS1FFT, imgS1FFTAmp
	  Image imgS1FFTAmpFiltered, filter, imgS2FFTAmp, imgS2FFTAmp0, imgS2FFTAmp_FFT, imgS2FFTAmp_FFTAmp, mask, ret
	  Image imgS1FFTAmpFiltered_Clip, imgS2FFTAmp_FFTAmp_Clip, imgS2FFTAmp0_Clip, imgS2FFTAmp_FFT_Ph
	  Image dst, dst0, projectionLines0, rotationAverage0, oneLine0, projectionLines, rotationAverage, oneLine
	  Image padded, paddedFFT, paddedFFTAmp, paddedFFTAmpHalf, padded1, paddedFFT1, paddedFFTAmp1, paddedFFTAmpHalf1
	  Image imgS2FFTShifted, imgS2FFTShiftedAmp, imgS2FFTShifted_FFTAmp, imgS2FFTShifted_FFTAmp2, dst2, dst2T, dst20, projectionLines2, rotationAverage2, oneLine2
	  Image oneLine2Minor, oneLine2Major, ret2, ellipseRadius, ellipseRadiusFindMax, oneLine3
	  Image ellipseProjectionLines, ellipseRotationAverage, ellipseOneLine, trajectory, trajectoryInv
	  Image img2, fftavg, fftavgAmp, bimg, tempfft, tempfftAmp
	  Image oneLineMinor, oneLineMajor, dst30
	  //Number size
	  
	  
     
     
	idcListenFFT(Object self)
    {
		num_chg = 0
		result("\n object 0x"+self.scriptObjectGetID().hex()+" created.")
    }
  
	~idcListenFFT(Object self)
    {
		result("\n object 0x"+self.scriptObjectGetID().hex()+" destroyed.")
    }


//correct anisotropic magnification distortion, see RefineAnisotropicScaleAligner
Image anisotropicScaling(Object self, Image img, Number aniso, Number theta, Number scale)
{
	Number dx = 0.0, dy = 0.0
	Number nx, ny
	getSize(img, nx, ny)
	result("\n scale "+scale+"\n aniso "+aniso+"\n thata "+theta+"\n nx*ny "+nx+"*"+ny)
	
	/*
	rotaton matrix that will rotate image counter-clock wise
	m11    m12
	m21    m22
	*/

	//3 steps:
	//first: rotate -theta
	//second: scale along x and y axes
	//thirs: rotate theta back to original orientation
	Number sx = 1./(scale*(1 + aniso*0.5))
	//Number sy = 1./(scale*(1 - aniso*0.5))
	Number sy = (1 + aniso*0.5)/scale //make the ellipse area constant for different aniso parameter
	sx = 1/sx
	sy = 1/sy

	Number cosine = cos(theta * pi() / 180.)
	Number sine = sin(theta * pi() / 180.)
	//result("\n sx, sy = "+sx+" "+sy)
	//result("\n cos, sin = "+cosine+" "+sine)

	Number m11 = sx*cosine*cosine+sy*sine*sine
	Number m22 = sx*sine*sine+sy*cosine*cosine
	Number m12 = (sx-sy)*cosine*sine
	Number m21 = m12
	//result("\n "+m11+" "+m12+"\n "+m21+" "+m22)

	Image correctedImage := RealImage("Corrected Image", 4, nx, ny)

	correctedImage = warp(img, \
					 m11*(icol-(nx/2+dx)) + m12*(irow-(ny/2+dy)) + (nx/2+dx), \
					 m21*(icol-(nx/2+dx)) + m22*(irow-(ny/2+dy)) + (ny/2+dy))


	//showImage(correctedImage)
	
	return correctedImage

}


//grid-boxing, 50% overlapping boxes
Image gridBoxing(Object self, Image img, Number boxsize)
	
	{
		
		//perform grid boxing with overlapping boxes
		result("\n Start grid boxing!")
		Number count = 0, nx, ny
		//Number apixX, apixY
		//apixX = imageGetDimensionScale( imgS1FFTAmp, 0 )
		//apixY = imageGetDimensionScale( imgS1FFTAmp, 1 )
		GetSize(img, nx, ny)
		
		//50% overlapping boxes
		Number cols = (floor(2*nx/boxsize/2))*2 - 1   //50% overlapping boxes
		Number rows = (floor(2*ny/boxsize/2))*2 - 1
		Number dx = boxsize/2
		Number dy = boxsize/2
		Number x0 = nx/2 - (floor(cols/2)+1)*boxsize/2
		Number y0 = ny/2 - (floor(rows/2)+1)*boxsize/2
		result("\n cols, rows = "+cols+", "+rows)
		Number mu, sigma, sigma0 = sqrt(variance(img))
		
		/*
		//75% overlapping boxes
		Number cols = (floor(2*nx/boxsize/2))*4 - 4   //75% overlapping boxes
		Number rows = (floor(2*ny/boxsize/2))*4 - 4
		Number dx = boxsize/4
		Number dy = boxsize/4
		Number x0 = nx/2 - (floor(cols/2)+1)*boxsize/4
		Number y0 = ny/2 - (floor(rows/2)+1)*boxsize/4
		result("\n "+cols+" "+rows+" "+x0+" "+y0)
		Number mu, sigma, sigma0 = sqrt(variance(img))
		*/
		
		
		//Number count = 0, nx, ny
		//GetSize(img, nx, ny)
		//Number mindim=min(nx, ny)
		//Number base2=trunc(log2(mindim))
		//Number boxsize = 2**(base2-1)
		result("\n ny*ny = "+nx+"*"+ny+"\n boxsize = "+boxsize)
		/*
		Number cols = floor(nx/boxsize)
		Number rows = floor(ny/boxsize)
		result("\n cols, rows = "+cols+", "+rows)

		Number x0 = floor((nx - cols*boxsize)/2), y0 = floor((ny - rows*boxsize)/2)
		result("\n x0, y0 = "+x0+", "+y0)

		Number dx = boxsize, dy = boxsize
		Number mu, sigma, sigma0 = sqrt(variance(img))
		//result("\n cols, rows = "+cols+", "+rows)
		//result("\n x0, y0 = "+x0+", "+y0)
		*/
		
		Image fftavg := complexImage("imgS1FFT_Complex", 8, boxsize, boxsize)
		Image fftavgAmp := realImage("imgS1FFT_Real", 4, boxsize, boxsize)
		Image tempfft := complexImage("tempfft_Complex", 8, boxsize, boxsize)
		Image tempfftAmp := realImage("tempfft_Real", 4, boxsize, boxsize)
		fftavg = 0
		fftavgAmp = 0
		
		Number rowStart = floor(rows/2-0), rowEnd = ceil(rows/2+1)
		Number colStart = floor(cols/2-1), colEnd = ceil(cols/2+1)
		//result("\n Grid boxing......rows:"+rowStart+"-"+rowEnd+" cols:"+colStart+"-"+colEnd)
		for (Number j=0; j<rows; j++)
		//for (Number j=rowStart; j<rowEnd; j++)
		{
			for (Number i=0; i<cols; i++)
			//for (Number i=colStart; i<colEnd; i++)
			{
				//result("\n "+i+" "+(x0+i*dx)+" "+j+" "+(y0+j*dy))
				Image bimg := img.Slice2(x0+i*dx, y0+j*dy, 0, 0, boxsize, 1, 1, boxsize, 1)
				//if (i ==0 && j == 0) showImage(bimg)
				mu = mean(bimg)
				sigma = sqrt(variance(bimg))
				//bimg = (bimg-mu)/sigma  //normalize
				
				//if (!(sigma0/3.0 < sigma < sigma0*3.0)) continue
				
				convertToFloat(bimg)
				tempfft = 0
				tempfftAmp = 0
				tempfft = RealFFT(bimg)/boxsize
				tempfftAmp = log10(modulus(tempfft))
				//if (i ==0 && j == 0) showImage(tempfftAmp)
				fftavgAmp = distance(tempfftAmp[icol, irow], fftavgAmp[icol, irow])
				
				//fftavg = distance(tempfft[icol, irow], fftavg[icol, irow])
				
				count += 1
			
			}
		}
		
		fftavgAmp /= sqrt(count)
		//fftavg/= sqrt(count)
		//ShowImage(fftavg)
		return fftavgAmp

	}


//pad 1D array
Image oversample1DArray(Object self, Number oversample, Number nx, Number nsamples, Image array)
	{
			Image padded := createFloatImage( "Padded_Projection_One_Line0", nsamples, 1)
			padded = 0
			//padded[0, 0, 1, floor((oversample/2-0.5)*nx)] = 0
			//padded[0, floor((oversample/2+0.5)*nx), 1, nsamples] = 0
			padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = array[icol, irow]
			Image paddedFFT := RealFFT(padded)
			Image paddedFFTAmp := (modulus(paddedFFT))
			Image paddedFFTAmp2 := paddedFFTAmp[0, nsamples/2, 1, nsamples]

			return paddedFFTAmp2 

	}
	
	
	
image GaussianConvolution(Object self, image sourceimg, number standarddev)
	{
	
		// Trap for zero (or negative) standard deviations
		
		if(standarddev<=0) return sourceimg
	
	
		// get the size of the source image. If it is not a power of 2 in dimension
		// warp it so that it is.

		number xsize, ysize, div2size, expandx, expandy, logsize
		getsize(sourceimg, xsize, ysize)
		expandx=xsize
		expandy=ysize


		// Check the x axis for power of 2 dimension - if it is not, round up to the next size
		// eg if it is 257 pixels round it up to 512.

		logsize=log2(xsize)
		if(mod(logsize,1)!=0) logsize=logsize-mod(logsize,1)+1
		expandx=2**logsize


		// Check the y axis for power of 2 dimension - if it is not, round up to the next size
		// eg if it is 257 pixels round it up to 512.

		logsize=log2(ysize)
		if(mod(logsize,1)!=0) logsize=logsize-mod(logsize,1)+1
		expandy=2**logsize


		// Use the Warp function to stretch the image to fit into the revised dimensions

		image warpimg=realimage("",4,expandx, expandy)
		warpimg=warp(sourceimg, icol*xsize/expandx, irow*ysize/expandy)


		// Create the gaussian kernel using the same dimensions as the expanded image

		image kernelimg:=realimage("",4,expandx,expandy)
		number xmidpoint=xsize/2
		number ymidpoint=ysize/2
		kernelimg=1/(2*pi()*standarddev**2)*exp(-1*(((icol-xmidpoint)**2+(irow-ymidpoint)**2)/(2*standarddev**2)))


		// Carry out the convolution in Fourier space

		compleximage fftkernelimg:=realFFT(kernelimg)
		compleximage FFTSource:=realfft(warpimg)
		compleximage FFTProduct:=FFTSource*fftkernelimg.modulus().sqrt()
		realimage invFFT:=realIFFT(FFTProduct)


		// Warp the convoluted image back to the original size

		image filter=realimage("",4,xsize, ysize)
		filter=warp(invFFT,icol/xsize*expandx,irow/ysize*expandy)
		return filter
	}
	
	
	
//calculate the mean defocus from the original single ring and calculate the max amount of shift for Fourier shifting
Image calculateShiftAndDefocus(Object self, Image imgS1FFTAmp, Number rmax, Number s2max, Number apix)
    {
		
		Number rangeDeg = 360, numWdg = 18
		Number apixX = apix, apixY = apix, numX, numY
		getSize(imgS1FFTAmp, numX, numY)
		result("\n In calculateShiftAndDefocus, apix "+apixX+"\n image size "+numX+"*"+numY) //here numX==numY
		//if (verbose == 2) ShowImage(imgS1FFTAmp)
		
		//add a mask in Fourier space to serve as a high pass filter to the input image 
		filter := RealImage( "Filter", 4, numX, numY )
		filter = tert(distance((icol-numX/2),(irow-numY/2))<(numX/2)*0.02, 0, 1) //set cutoff fixed?
		imgS1FFTAmpFiltered = RealIFFT( RealFFT(imgS1FFTAmp) * filter )
		imageSetName(imgS1FFTAmpFiltered, "imgS1_FFT_Amplitude_Filtered")
		if (clip && verbose > 1) 
			{
				ShowImage(imgS1FFTAmpFiltered) //2D-S1 power spectra
				Number startX = round((boxsize-clip)/2), startY = round((boxsize-clip)/2)
				imgS1FFTAmpFiltered_Clip := RealImage( "imgS1_FFT_Amplitude_Filtered_Clip", 4, clip, clip )
				imgS1FFTAmpFiltered_Clip = imgS1FFTAmpFiltered.Slice2(startX, startY, 0, 0, clip, 1, 1, clip, 1)
				ShowImage(imgS1FFTAmpFiltered_Clip)
				
			}
		
		
		result("\n Convert S1-->S2")
		imgS2FFTAmp := RealImage( "imgS2_FFT_Amplitude", 4, numX, numY )
		imgS2FFTAmp0 := RealImage( "imgS2_FFT_Amplitude", 4, numX, numY )
		imgS2FFTAmp = warp( imgS1FFTAmpFiltered, \
		sqrt(distance((icol-numX/2),(irow-numY/2))/(numX/2))*rmax*cos(atan2(irow-numY/2, icol-numX/2))+numX/2,\
		sqrt(distance((icol-numX/2),(irow-numY/2))/(numY/2))*rmax*sin(atan2(irow-numY/2, icol-numX/2))+numY/2 )
		
		imgS2FFTAmp0 = imgS2FFTAmp.imageClone()
		if (clip && verbose > 1) 
			 {
				ShowImage(imgS2FFTAmp0)  //2D-S2 power spectra
				Number startX = round((boxsize-clip)/2), startY = round((boxsize-clip)/2)
				imgS2FFTAmp0_Clip := RealImage( "imgS2_FFT_Amplitude_Clip", 4, clip, clip )
				imgS2FFTAmp0_Clip = imgS2FFTAmp0.Slice2(startX, startY, 0, 0, clip, 1, 1, clip, 1)
				ShowImage(imgS2FFTAmp0_Clip)
				
			 }
		//showimage(imgS2FFTAmp)
		
		/*
		Image gaussblur := self.GaussianConvolution(imgS2FFTAmp,5)
		Image imgS2FFTAmp_FFT_gauss := ComplexImage("imgS2_FFT_Amplitude_FFT_gauss", 8, numX, numY)
		imgS2FFTAmp_FFT_gauss = RealFFT(gaussblur)
		ShowImage(imgS2FFTAmp_FFT_gauss)
		*/
		Number halfMinor = min(numX, numY)/2
		Number nsamples = halfMinor*(oversample), nx = halfMinor
		//Number peak = 3*2*s2max*10000*wavelen
		//result("\n peak = "+peak)
		filter = tert(distance((icol-numX/2),(irow-numY/2))<(numX/2)*0.05*(1/apixX), 0, 1)
		imgS2FFTAmp := RealIFFT( RealFFT(imgS2FFTAmp) * filter )
		
		
		
		Number thresh
		thresh = max(min(imgS2FFTAmp), mean(imgS2FFTAmp)-3*sqrt(variance(imgS2FFTAmp)))
		imgS2FFTAmp = tert(imgS2FFTAmp[icol, irow]<thresh, thresh, imgS2FFTAmp[icol, irow])
		//imageSetName(imgS2FFTAmp, "imgS2_FFT_Amplitude")
		
		
		imgS2FFTAmp_FFT := ComplexImage("imgS2_FFT_Amplitude_FFT", 8, numX, numY)
		imgS2FFTAmp_FFT = RealFFT(imgS2FFTAmp)
		


		imgS2FFTAmp_FFTAmp := RealImage( "imgS2_FFT_Amplitude_FFT_Amplitude", 4, numX, numY )
		imgS2FFTAmp_FFTAmp = log10(modulus(imgS2FFTAmp_FFT)) //will be transferred to the main function and used later
		
		//imgS2FFTAmp_FFT_Ph := RealImage( "imgS2_FFT_Amplitude_FFT_Phase", 4, numX, numY )
		//imgS2FFTAmp_FFT_Ph = Phase(imgS2FFTAmp_FFT)
		//ShowImage(imgS2FFTAmp_FFT_Ph)
		//ShowImage(imgS2FFTAmp_FFTAmp)
		//imgS2FFTAmp_FFTAmp = (modulus(imgS2FFTAmp_FFT))
		
		//imageSetName(imgS2FFTAmp_FFTAmp, imageName)
		
		if (clip && verbose > 1)
		{
			ShowImage(imgS2FFTAmp_FFTAmp)			
			Number startX = round((boxsize-clip)/2), startY = round((boxsize-clip)/2)
			imgS2FFTAmp_FFTAmp_Clip := RealImage( "imgS2_FFT_Amplitude_FFT_Amplitude_Clip", 4, clip, clip )
			imgS2FFTAmp_FFTAmp_Clip = imgS2FFTAmp_FFTAmp.Slice2(startX, startY, 0, 0, clip, 1, 1, clip, 1)
			ShowImage(imgS2FFTAmp_FFTAmp_Clip)
		}

		//add a mask to serve as high pass filter, Fourier shifted image may not need this mask
		mask := RealImage( "mask", 4, numX, numY )
		mask = tert(distance((icol-numX/2),(irow-numY/2))<(numX/2)*0.02*(1/apixX), 0, 1) //set cutoff fixed?
		imgS2FFTAmp = imgS2FFTAmp*mask
    
		//find the center of the image
		Number centerX=numX/2, centerY = numY/2
		//determine the half of the smallest dimension
		//Number halfMinor = min(numX, numY)/2
		
		Number dr = 2*pi() / rangeDeg //if rangeDeg = 360, dr = 2*pi()/rangeDeg
		Number polarSamples = rangeDeg
		//convert the image to polar coordinates dst
		//they will be transferred to the main function and used later
		//dst is the polar coordinates of imgS2FFTAmp_FFTAmp
		//dst0 is the polar coordinates of imgS2FFTAmp
		dst := createFloatImage( "dst_polar_coordinates_imgS2FFT_FFTAmp", halfMinor, polarSamples )
		dst0 := createFloatImage( "dst_polar_coordinates_imgS2FFTAmp", halfMinor, polarSamples )
		dst = warp( imgS2FFTAmp_FFTAmp, icol*sin(irow*dr) + centerX, icol*cos(irow*dr) + centerY )
		dst0 = warp( imgS2FFTAmp, icol*sin(irow*dr) + centerX, icol*cos(irow*dr) + centerY )
		//ShowImage(dst)
		if (verbose > 1)
		{
			Image dstT := matrixTranspose(dst)
			dstT.setName("dst_polar_coordinates_imgS2FFT_FFTAmp_transposed")
			FlipVertical(dstT)
			ShowImage(dstT)
		}
		
		/*
		//projectionLines: roational average of each wedge, # of lines=numWdg
		projectionLines := createFloatImage( "Projection_Lines", halfMinor, numWdg )
		Image oneLine := createFloatImage("Projection_One_Line", halfMinor, 1)
		//rotationAverage: store the lines of each wedge
		Image rotationAverage := createFloatImage( "Rotational_Average", halfMinor, rangeDeg/numWdg )
		projectionLines = 0
		*/
		//Number oversample = 100 //should pad the 1D lines extracted from dst0 (imgS2FFTAmp)
		//Number nsamples = halfMinor*(oversample), nx = halfMinor
		Number shift = 0, defocus = 0, peakPosition = 0
		
		Image projectionLines0 := createFloatImage( "Projection_Lines0", nsamples/2, numWdg )
		Image oneLine0 := createFloatImage("Projection_One_Line0", halfMinor, 1)
		Image rotationAverage0 := createFloatImage( "Rotational_Average0", halfMinor, rangeDeg/numWdg )
		projectionLines0 = 0
		
		Number r_rgb = 0, g_rgb = 1, b_rgb = 0
		for (Number curWdg=0; curWdg<numWdg; curWdg++)
			{
			
			r_rgb = (Mod(curWdg,4)==0)?(Abs(r_rgb-1)):r_rgb
			g_rgb = (Mod(curWdg,2)==0)?(Abs(g_rgb-1)):g_rgb
			b_rgb = (Mod(curWdg,1)==0)?(Abs(b_rgb-1)):b_rgb
			if (r_rgb + g_rgb + b_rgb == 3)
				{
				r_rgb = 0
				g_rgb = 0
				b_rgb = 0
				}

			Number angle = curWdg*(RangeDeg/numWdg)

			//rotationAverage = 0
			rotationAverage0 = 0
			
			//rotationAverage = dst[curWdg*(RangeDeg/numWdg), 0, (curWdg+1)*(RangeDeg/numWdg), halfMinor]
			rotationAverage0 = dst0[curWdg*(RangeDeg/numWdg), 0, (curWdg+1)*(RangeDeg/numWdg), halfMinor]
			
			//oneLine = 0
			oneLine0 = 0
			
			//oneLine[icol, irow] += rotationAverage
			oneLine0[icol, irow] += rotationAverage0
			oneLine0 /= (RangeDeg/numWdg) //normalize
			
			//projectionLines[curWdg, 0, curWdg+1, halfMinor] = oneLine
			//projectionLines[curWdg, 0, curWdg+1, halfMinor] /= (RangeDeg/numWdg) //normalize
			
			//Number maxX0, maxY0
			//max(oneLine, maxX0, maxY0)
			//result("\n Current degree "+angle+" Find peak position in 1D: " + "x0 = "+maxX0) 
			
			Image paddedFFTAmp2 := self.oversample1DArray(oversample, nx, nsamples, oneLine0) 
			//projectionLines0[curWdg, 0, curWdg+1, nsamples/2] = paddedFFTAmp2 
			
			Number dfX0, dfY0
			max(paddedFFTAmp2, dfX0, dfY0)
			//max(paddedFFTAmp, dfX0, dfY0)
			//result("\n Current degree "+angle+" Find peak position in 1D: " + "x0 = "+dfX0) 
			Number df = dfX0 /((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
			
			defocus += df
			shift += (nsamples/2 - dfX0)/(oversample/2)
			peakPosition += dfX0/(oversample/2)
			
			//result("\n defocus: "+df+"um"+" Find peak position in padded 1D: "+"x0 = "+dfX0+"  maxdShift = "+(nsamples/2 - dfX0)/(oversample/2))
			/*
			setDisplayType(projectionLines, 4)
			linePlotImageDisplay lpid = projectionLines.imageGetImageDisplay(0)
			lpid.linePlotImageDisplaySetDoAutoSurvey(0,0)
			lpid.linePlotImageDisplaySetGridOn(0)
			lpid.linePlotImageDisplaySetSliceDrawingStyle(curWdg, 1)
			lpid.lineplotimagedisplaysetcontrastlimits(min(projectionLines)-0.5, max(projectionLines)+0.5)
			
			setDisplayType(projectionLines0, 4)
			linePlotImageDisplay lpid0 = projectionLines0.imageGetImageDisplay(0)
			lpid0.linePlotImageDisplaySetDoAutoSurvey(0,0)
			lpid0.linePlotImageDisplaySetGridOn(0)
			lpid0.linePlotImageDisplaySetSliceDrawingStyle(curWdg, 1)
			lpid0.lineplotimagedisplaysetcontrastlimits(min(projectionLines0)-0.5, max(projectionLines0)+0.5)
			*/
		}
		
		shift /= numWdg
		defocus /= numWdg
		peakPosition /= numWdg
		result("\n mean shift = "+shift+" ; mean defocus = "+defocus+"um")
		
		Image ret := RealImage( "recommendShift&meanDefocus", 4, 3, 1 )
		ret[0, 0, 1, 1] = shift
		ret[0, 1, 1, 2] = defocus
		ret[0, 2, 1, 3] = peakPosition
		
		//displayAt(projectionLines0, 139, 47)
		//setWindowSize(projectionLines0, 800, 500)
		
		//displayAt(projectionLines, 139, 47)
		//setWindowSize(projectionLines, 800, 500)
		result("\n Finish calculatedShiftAndDefocus")
		
		return ret

	}




//pad 1D array, multiply exp((-1)*j*(2*pi()*shift*oversample/2) to perform 1D Fourier shifting, then FFT
Image oversample1DArrayShifted(Object self, Number oversample, Number nx, Number nsamples, Image array, Number shift)
	{

			Image padded := createFloatImage( "Padded_Projection_One_Line0", nsamples, 1 )
			padded = 0
			//padded[0, 0, 1, floor((oversample/2-0.5)*nx)] = 0
			//padded[0, floor((oversample/2+0.5)*nx), 1, nsamples] = 0
			padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = array[icol, irow]
			
			Image paddedShifted := ComplexImage( "Padded_Projection_One_Line0_Shifted", 8, nsamples, 1 )
			ComplexNumber j = complex(0, 1) 
			paddedShifted = padded[icol, irow] * exp((-1)*j*(2*pi()*shift*oversample/2)*(icol-nsamples/2) / (nsamples) )
			Image paddedFFT := FFT(paddedShifted)
			//Image paddedFFT := RealFFT(padded)
			Image paddedFFTAmp := (modulus(paddedFFT))
			Image paddedFFTAmp2 := paddedFFTAmp[0, nsamples/2, 1, nsamples] 
			
			return paddedFFTAmp2


}


Number refinePhi(Object self, Number oversample, Number nx, Number nsamples, Image dst30, Number shift, Image mask1D, Number start, Number end)
{	
			
		result("\n startPhi = "+start+"  endPhi = "+end)
		Number dst30X, dst30Y
		getSize(dst30, dst30X, dst30Y)
		result("\n dst30 image size "+dst30X+"*"+dst30Y)
		//In dst30, each line means 1/fineFactor degree in 360 degree
		Number maxCurWdg = 0
		Number maxX2 = 0, maxY2 = 0, maxX = 0
		Number step = 1
		
		for (Number curWdg=start*fineFactor3; curWdg<end*fineFactor3; curWdg+=step)
		{
			
				Image oneLine0 := dst30[curWdg, 0, (curWdg+1), nx]
				Image paddedFFTAmp2 := self.oversample1DArrayShifted(oversample, nx, nsamples, oneLine0, shift) 
				paddedFFTAmp2 = paddedFFTAmp2*mask1D
				max(paddedFFTAmp2, maxX2, maxY2)
				if (maxX2 > maxX)
				{
					maxX = maxX2
					maxCurWdg = curWdg
				}
				
		}
		
		
		Number refinedPhi = maxCurWdg/fineFactor3
		result("\n refinedPhi = "+refinedPhi)
		
		return refinedPhi


}

//calculate the major/minor direction of the ellipse
Image calculateMinorMajorDirection(Object self, Image imgS2FFTAmp, Number shift, Number peakPosition, Number s2max, Number apix)
	{
		
		Number rangeDeg = 360, numWdg = numWdg1
		Number apixX = apix, apixY = apix, numX, numY
		getSize(imgS2FFTAmp, numX, numY)
		result("\n In calculateMinorMajorDirection, apix "+apixX+" "+apixY+"\n image size "+numX+"*"+numY) //here numX==numY
		
		imgS2FFTShifted := ComplexImage("imgS2_FFT_Shifted", 8, numX, numY )
		ComplexNumber j = complex(0, 1) 
		result("\n Real shift = "+shift)
		imgS2FFTShifted = imgS2FFTAmp[icol, irow] * exp(j * 2* pi()* shift * distance((icol-numX/2),(irow-numY/2)) / (numX))
		
		imgS2FFTShifted_FFTAmp := RealImage("imgS2_FFT_Shifted_FFT_Amplitude", 4, numX, numY )
		imgS2FFTShifted_FFTAmp = log10(modulus(FFT(imgS2FFTShifted)))
		
		
		//How large the mask is? Do padding to remove the mask?
		mask := RealImage( "mask", 4, numX, numY )
		//mask = tert(distance((icol-numX/2),(irow-numY/2)) < numX/2*0.2*(1/apixX), 0, 1) //set cutoff fixed?
		mask = tert(distance((icol-numX/2),(irow-numY/2)) < peakPosition+shift-15 || distance((icol-numX/2),(irow-numY/2)) > peakPosition+shift+15, 0, 1) //we can calculate the mask size
		//mask = tert( distance((icol-numX/2),(irow-numY/2)) > peakPosition+10 && distance((icol-numX/2),(irow-numY/2)) < peakPosition+40, 1, 0) //set cutoff fixed?
		imgS2FFTShifted_FFTAmp = imgS2FFTShifted_FFTAmp*mask
		if (verbose > 1) ShowImage(imgS2FFTShifted_FFTAmp) //Fourier shifted 2D-S2 power spectra
		//ShowImage(imgS2FFTShifted_FFTAmp)
		
		//find the center of the image
		Number centerX = numX/2, centerY = numY/2
		//determine the half of the smallest dimension
		Number halfMinor = min(numX, numY)/2
		Number polarSamples2 = rangeDeg * fineFactor2
		Number polarSamples3 = rangeDeg * fineFactor3
		Number dr = 2*pi() / rangeDeg
		Number dr2 = 2*pi() / polarSamples2 //if rangeDeg = 360, dr = 2*pi()/rangeDeg
		Number dr3 = 2*pi() / polarSamples3
		
		
		//convert the image to polar coordinates
		dst2 := createFloatImage( "dst2_Polar_Coordinates_imgS2FFTShifted_FFTAmp", halfMinor, polarSamples2 )
		dst2 = warp( imgS2FFTShifted_FFTAmp, icol*sin(irow*dr2) + centerX, icol*cos(irow*dr2) + centerY )
		if (verbose > 1)
		{
			dst2T := matrixTranspose(dst2)
			dst2T.setName("dst2_Polar_Coordinates_imgS2FFTShifted_FFTAmp_transposed")
			FlipVertical(dst2T)
			ShowImage(dst2T) //polar cooredinates of Fourier shifted 2D-S2 power spectra
		}
		
		dst20 := createFloatImage( "dst20_Polar_Coordinates", halfMinor, polarSamples2 )
		dst30 := createFloatImage( "dst30_Polar_Coordinates", halfMinor, polarSamples3 ) //finer sampling in polar coordinate, will use dst30 in calculateMinorandMajorRadius
		//mask = tert(distance((icol-numX/2),(irow-numY/2))<(numX/2)*0.04*(1/apixX), 0, 1)
		//imgS2FFTAmp = imgS2FFTAmp*mask
		dst20 = warp( imgS2FFTAmp, icol*sin(irow*dr2) + centerX, icol*cos(irow*dr2) + centerY )
		dst30 = warp( imgS2FFTAmp, icol*sin(irow*dr3) + centerX, icol*cos(irow*dr3) + centerY )
		//dst20 = warp( imgS2FFTShiftedAmp, icol*sin(irow*dr2) + centerX, icol*cos(irow*dr2) + centerY )
		//ShowImage(dst20)

		
		Number x0 = peakPosition
		result("\n Initial detection of peak position in polar coordinate: "+x0)
		Number halfMinorNew = (x0-halfMinor/2)>0?(halfMinor-x0)*2:x0*2
		Number offset = halfMinor - halfMinorNew
		Number step = 1

		//projectionLines: roational average of each wedge, # of lines=numWdg
		Image projectionLines2 := createFloatImage( "Projection_Lines2", halfMinor, numWdg )
		Image oneLine2 := createFloatImage("Projection_One_Line2", halfMinor, 1)
		projectionLines2 = 0
		//rotationAverage: store the lines of each wedge
		//Image rotationAverage2 := createFloatImage( "Rotational_Average2", halfMinor, polarSamples2/numWdg )
		Image ellipseRadius := createFloatImage("Ellipse_Radius", numWdg/step, 1)
		
		Number peakPositionOversample = (peakPosition+shift) * (oversample/2)
		//Number peakPositionOversample = (peakPosition) * (oversample/2)
		result("\n peakPositionOversample = " +peakPositionOversample)
		result("\n shift*oversample/2 = " +shift*oversample/2)
		Number nsamples = halfMinor*(oversample), nx = halfMinor
		Image mask1D := RealImage( "mask", 4, nsamples/2, 1 )
		//mask1D = tert(icol<(peakPosition+shift+15)*(oversample/2) && icol>(peakPosition+shift-15)*(oversample/2), 1, 0)
		mask1D = tert(icol<peakPositionOversample*0.9 || icol>peakPositionOversample*1.1, 0, 1)
		result("\n "+peakPositionOversample*0.9+" "+peakPositionOversample*1.1)
		
		projectionLines0 := createFloatImage( "Projection_Lines0", nsamples/2, numWdg )
		Image oneLine0 := createFloatImage("Projection_One_Line0", halfMinor, 1)
		//Image rotationAverage0 := createFloatImage( "Rotational_Average0", halfMinor, samples/numWdg )
		projectionLines0 = 0
		
		Number r_rgb = 0, g_rgb = 1, b_rgb = 0
		for (Number curWdg=0; curWdg<numWdg; curWdg+=step)
		{
			
			r_rgb = (Mod(curWdg,4)==0)?(Abs(r_rgb-1)):r_rgb
			g_rgb = (Mod(curWdg,2)==0)?(Abs(g_rgb-1)):g_rgb
			b_rgb = (Mod(curWdg,1)==0)?(Abs(b_rgb-1)):b_rgb
			if (r_rgb + g_rgb + b_rgb == 3)
				{
				r_rgb = 0
				g_rgb = 0
				b_rgb = 0
				}

			//Number angle = curWdg*(rangeDeg/numWdg)
			
			//Image rotationAverage2 := dst2[(curWdg)*(polarSamples2/numWdg), 0, (curWdg+1)*(polarSamples2/numWdg), halfMinor]
			Image rotationAverage0 := dst20[curWdg*(polarSamples2/numWdg), 0, (curWdg+1)*(polarSamples2/numWdg), halfMinor]
			//Image rotationAverage0 := dst20[curWdg*(polarSamples2/numWdg), 0, (curWdg)*(polarSamples2/numWdg)+1, halfMinor]
			
			//oneLine2 = 0
			oneLine0 = 0
			
			//oneLine2[icol, irow] += rotationAverage2
			oneLine0[icol, irow] += rotationAverage0
			oneLine0 /= (polarSamples2/numWdg) //normalizes
			
			//projectionLines2[curWdg, 0, curWdg+1, halfMinor] = oneLine2
			//projectionLines2[curWdg, 0, curWdg+1, halfMinor] /= (polarSamples2/numWdg) //normalize
			
			Image paddedFFTAmp2 := self.oversample1DArrayShifted(oversample, nx, nsamples, oneLine0, shift) 
			paddedFFTAmp2 = paddedFFTAmp2*mask1D
			//projectionLines0[curWdg, 0, curWdg+1, nsamples/2] = paddedFFTAmp2
	
			//oneLine2 = dst2[icol, curWdg]
			Number maxX2, maxY2
			//max(oneLine2, maxX2, maxY2)
			max(paddedFFTAmp2, maxX2, maxY2)
			maxX2 /= (oversample/2)
			Number cur = curWdg/step
			ellipseRadius[0, cur, 1, cur+1] = maxX2 - shift //remove the amount of Fourier shift
			//result("\n "+curWdg+" "+maxX2)
			
			/*
			setDisplayType(projectionLines2, 4)
			projectionLines2.ImageSetDimensionUnitString(0, "Radius (pixels)" )
			projectionLines2.ImageSetIntensityUnitString("Intensity")
			linePlotImageDisplay lpid2 = projectionLines2.imageGetImageDisplay(0)
			lpid2.linePlotImageDisplaySetDoAutoSurvey(0,0)
			lpid2.linePlotImageDisplaySetGridOn(1)
			lpid2.linePlotImageDisplaySetSliceDrawingStyle(curWdg, 1)
			//lpid2.lineplotimagedisplaysetdisplayedchannels(offset, halfMinor)
			lpid2.lineplotimagedisplaysetcontrastlimits(min(projectionLines2)-0.5, max(projectionLines2)+0.5)
			
			
			setDisplayType(projectionLines0, 4)
			projectionLines0.ImageSetDimensionUnitString(0, "Radius (pixels)" )
			projectionLines0.ImageSetIntensityUnitString("Intensity")
			linePlotImageDisplay lpid0 = projectionLines0.imageGetImageDisplay(0)
			lpid0.linePlotImageDisplaySetDoAutoSurvey(0,0)
			lpid0.linePlotImageDisplaySetGridOn(1)
			lpid0.linePlotImageDisplaySetSliceDrawingStyle(curWdg, 1)
			//lpid0.lineplotimagedisplaysetdisplayedchannels(offset, halfMinor)
			lpid0.lineplotimagedisplaysetcontrastlimits(min(projectionLines0)-0.5, max(projectionLines0)+0.5)
			*/
			
			
			
		}
		
		//fix the extremely small values due to the cross in power spectra
		Number eRadiusX, eRadiusY
		Number eRadiusMean = mean(ellipseRadius)
		eRadiusMean = mean(ellipseRadius)
		getSize(ellipseRadius, eRadiusX, eRadiusY)
		for (Number i=1; i<eRadiusX-1; i++)
		{
			if (sum(ellipseRadius[0, i, 1, i+1]) < (eRadiusMean/2.0))
			{
				Number pre = sum(ellipseRadius[0, i-1, 1, i])
				Number post = sum(ellipseRadius[0, i+1, 1, i+2])
				ellipseRadius[0, i, 1, i+1] = (pre+post)/2.0
			}
		}

		
		
		
		
		setDisplayType(ellipseRadius, 4)
		ImageSetDimensionScale(ellipseRadius, 0,  rangeDeg/numWdg)
		//convert Y-axis label from pixel to nm
		Number factor = 1000*(oversample/2)/((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		//factor = 1 
		ImageSetIntensityScale(ellipseRadius, factor)
		//ImageSetIntensityOrigin(ellipseRadius, shift)
		
		ellipseRadius.ImageSetDimensionUnitString(0, "Angle (degree)" )
		//ellipseRadius.ImageSetIntensityUnitString("Radius (pixels)" )
		ellipseRadius.ImageSetIntensityUnitString("Defocus (nm)" )
		linePlotImageDisplay lpid3 = ellipseRadius.imageGetImageDisplay(0)
		lpid3.linePlotImageDisplaySetDoAutoSurvey(0,0)
		Number minLimit = min(ellipseRadius), maxLimit = max(ellipseRadius), deltaLimit = maxLimit-minLimit
		lpid3.lineplotimagedisplaysetcontrastlimits(minLimit - deltaLimit*0.25, maxLimit + deltaLimit*0.25)
		//lpid3.lineplotimagedisplaysetfontsize(15)
		
		//calculate the major and minor axis of ellipse
		numWdg = numWdg/step
		ComplexImage ellipseRadiusFFT := ComplexImage("Ellipse_Radius_FFT", 8, numWdg, 1)
		ellipseRadiusFFT = RealFFT(ellipseRadius-mean(ellipseRadius))
		result("\n mean of ellipseRadius = "+ mean(ellipseRadius)*factor)
		
		setDisplayType(ellipseRadiusFFT, 4)
		//ImageSetDimensionScale(ellipseRadiusFFT, 0,  (2*pi()/numWdg))
		//ImageSetDimensionOrigin(ellipseRadiusFFT, 0, (-1)*(2*pi()/numWdg)*(numWdg/2))
		ImageSetDimensionOrigin(ellipseRadiusFFT, 0, (-1)*(numWdg/2))
		//ellipseRadiusFFT.ImageSetDimensionUnitString(0, "Frequency (radian/sample)" )
		ellipseRadiusFFT.ImageSetDimensionUnitString(0, "Angular Frequency (Hz)" )
		//ellipseRadiusFFT.ImageSetIntensityUnitString( "Intensity" )
		ellipseRadiusFFT.ImageSetIntensityUnitString( "Amplitude" )
		linePlotImageDisplay lpid4 = ellipseRadiusFFT.imageGetImageDisplay(0)
		lpid4.linePlotImageDisplaySetDoAutoSurvey(0,0)
		//lpid4.lineplotimagedisplaysetfontsize(15)
		
		Number len = 10 //first 10 pixels
		ComplexImage ellipseRadiusFFTSub := ComplexImage("Ellipse_Radius_FFT_Sub", 8, len, 1) //Only 0~10 pixels of ellipseRadiusFFT
		ellipseRadiusFFTSub[0, 0, 1, 10] = ellipseRadiusFFT[0, numWdg/2, 1, numWdg/2+10]
		setDisplayType(ellipseRadiusFFTSub, 4)
		ellipseRadiusFFTSub.ImageSetDimensionUnitString(0, "Angular Frequency (Hz)" )
		//ellipseRadiusFFTSub.ImageSetIntensityUnitString( "Intensity" )
		ellipseRadiusFFTSub.ImageSetIntensityUnitString( "Amplitude" )
		linePlotImageDisplay lpid5 = ellipseRadiusFFTSub.imageGetImageDisplay(0)
		lpid5.linePlotImageDisplaySetDoAutoSurvey(0,0)
		//lpid5.lineplotimagedisplaysetfontsize(15)
		
		RealImage ellipseRadiusFFTSubMod := RealImage("Ellipse_Radius_FFT_Sub_Mod", 4, len, 1)
		ellipseRadiusFFTSubMod = modulus(ellipseRadiusFFTSub)
		Number high = max(ellipseRadiusFFTSubMod)
		Number low = 0
		//result("\n "+high+" "+low)
		//showImage(ellipseRadiusFFTSubMod)
		
		//make the x-axis label (e.g. 2Hz) below the bar positions
		//ComplexImage ellipseRadiusFFTSubModified := ComplexImage("Ellipse_Radius_FFT_Sub_Modified", 8, len*10, 1)
		RealImage ellipseRadiusFFTSubModified := REalImage("Ellipse_Radius_FFT_Sub_Modified", 4, len*10, 1)
		ellipseRadiusFFTSubModified = 0
		ellipseRadiusFFTSubModified[0, 5, 1, 15] = sum(ellipseRadiusFFTSubMod[0, 1, 1, 2])
		ellipseRadiusFFTSubModified[0, 15, 1, 25] = sum(ellipseRadiusFFTSubMod[0, 2, 1, 3])
		ellipseRadiusFFTSubModified[0, 25, 1, 35] = sum(ellipseRadiusFFTSubMod[0, 3, 1, 4])
		ellipseRadiusFFTSubModified[0, 35, 1, 45] = sum(ellipseRadiusFFTSubMod[0, 4, 1, 5])
		ellipseRadiusFFTSubModified[0, 45, 1, 55] = sum(ellipseRadiusFFTSubMod[0, 5, 1, 6])
		ellipseRadiusFFTSubModified[0, 55, 1, 65] = sum(ellipseRadiusFFTSubMod[0, 6, 1, 7])
		ellipseRadiusFFTSubModified[0, 65, 1, 75] = sum(ellipseRadiusFFTSubMod[0, 7, 1, 8])
		ellipseRadiusFFTSubModified[0, 75, 1, 85] = sum(ellipseRadiusFFTSubMod[0, 8, 1, 9])
		ellipseRadiusFFTSubModified[0, 85, 1, 95] = sum(ellipseRadiusFFTSubMod[0, 9, 1, 10])
		
		setDisplayType(ellipseRadiusFFTSubModified, 4)
		ImageSetDimensionScale(ellipseRadiusFFTSubModified, 0,  1/10)
		ellipseRadiusFFTSubModified.ImageSetDimensionUnitString(0, "Angular Frequency (Hz)" )
		//ellipseRadiusFFTSubModified.ImageSetIntensityUnitString( "Intensity" )
		ellipseRadiusFFTSubModified.ImageSetIntensityUnitString( "Amplitude" )
		linePlotImageDisplay lpid6 = ellipseRadiusFFTSubModified.imageGetImageDisplay(0)
		lpid6.linePlotImageDisplaySetDoAutoSurvey(0,0)
		//lpid6.lineplotimagedisplaysetfontsize(15)
		lpid6.lineplotimagedisplaySetContrastLimits(0, 10*round(high/10.0+0.5))
		//showimage(ellipseRadiusFFTSubModified)
		
		
		
		
		ComplexImage fitEllipseRadiusFFT := ComplexImage("Fit_Ellipse_Radius_FFT", 8, numWdg, 1)
		fitEllipseRadiusFFT = 0
		fitEllipseRadiusFFT[0, numWdg/2+2, 1, numWdg/2+3] = ellipseRadiusFFT[0, numWdg/2+2, 1, numWdg/2+3]
		fitEllipseRadiusFFT[0, numWdg/2-2, 1, numWdg/2-1] = ellipseRadiusFFT[0, numWdg/2-2, 1, numWdg/2-1]
		//showImage(fitEllipseRadiusFFT)
		
		image fitEllipseRadius := createFloatImage("Fit_Ellipse_Radius", numWdg, 1)
		fitEllipseRadius = 0
		fitEllipseRadius = RealIFFT(fitEllipseRadiusFFT)
		fitEllipseRadius += mean(ellipseRadius)
		//showImage(fitEllipseRadius)
		
		//lpid3.lineplotimagedisplaysetlegendshown(1)
		lpid3.imagedisplayaddimage(fitEllipseRadius, "Fitted Ellipse")
		object sliceid0=lpid3.imagedisplaygetsliceidbyindex(0)
		lpid3.imagedisplaysetslicelabelbyid(sliceid0, "Ellipse Radius")
		object sliceid1=lpid3.imagedisplaygetsliceidbyindex(1)
		if (GMS>1) lpid3.lineplotimagedisplaysetsliceLineThickness(sliceid1, 2)
		
		ComplexImage c = ellipseRadiusFFT[0, numWdg/2+2, 1, numWdg/2+3]
		ConvertToComplex(c)
		//Number phi2 = phase(c)
		Number re = sum(real(c)), im = sum(imaginary(c))
		ComplexNumber z = complex(re, im)
		Number phi = phase(z)*180/pi()  //in degree
		Number amp = abs(z)   //absolute value of the modulus of a complex number
		result("\n complexNumber "+re+"+"+im+"i "+"\n phase = "+phi+" degree"+" amp = "+amp)
		//calculate the length of major and minor axis from FFT directly
		Number calcMajor = mean(fitEllipseRadius)+amp/(numWdg/2)
		Number calcMinor = mean(fitEllipseRadius)-amp/(numWdg/2)
		result("\n calcMajor = "+calcMAjor+"  calcMinor = "+calcMinor)
		result("\n Max in Sine = "+max(fitEllipseRadius)+"  Min in Sine = "+min(fitEllipseRadius))
		
		//phi = mod(phi+360, 180)
		phi = phi/2  // 2-fold astigmatism
		phi = phi>0?phi:(phi+180)
		
		//Number start = floor(phi-(rangeDeg/numWdg)/2), end = ceil(phi+(rangeDeg/numWdg)/2)
		//result("\n startPhi = "+start+"  endPhi = "+end)
		//Number refinedPhi = self.refinePhi(oversample, nx, nsamples, dst30, shift, mask1D, start, end)
		//phi = refinedPhi
		
		
		Number majorDirection = phi
		//calculate the radius of phi direction
		Number curCol = floor(majorDirection*(polarSamples2/rangeDeg))
		Number majorRadius, minorRadius, temp1, temp2
		oneLine2Major := dst2[curCol, 0, curCol+1, halfMinor]
		max(oneLine2Major, majorRadius, temp1)
		
		//Number majorDirection = phi+90
		Number minorDirection = phi-90 >0 ? phi-90:phi+90
		curCol = floor(minorDirection*(polarSamples2/rangeDeg))
		//oneLine2Minor := dst20[curCol, 0, curCol+1, halfMinor]
		oneLine2Minor := dst2[curCol, 0, curCol+1, halfMinor]
		//showimage(oneLine2Minor)
		max(oneLine2Minor, minorRadius, temp2)
	

		//measure the major and minor axis length from the fitted curve
		Number minorRadius1, majorRadius1
		temp1 = majorDirection/360*numWdg
		temp2 = minorDirection/360*numWdg
		majorRadius1 = mean(fitEllipseRadius[0, temp1, 1, temp1+1])
		minorRadius1 = mean(fitEllipseRadius[0, temp2, 1, temp2+1])
		
		
		//majorRadius1 = calcMajor
		//minorRadius1 = calcMinor
		
		Number deltaDefocus, defocusMajor, defocusMinor
		deltaDefocus = (majorRadius1-minorRadius1)*factor
		defocusMajor = majorRadius1 * factor
		defocusMinor = minorRadius1 * factor
		
		
		
		result("\n majorDirection, minorDirection = "+majorDirection+" "+minorDirection)
		//result("\n majorRadius, minorRadius, diff = "+majorRadius+" "+minorRadius+" "+(majorRadius-minorRadius))
		result("\n From fitted curve, majorRadius, minorRadius, diff = "+majorRadius1+" "+minorRadius1+" "+(majorRadius1-minorRadius1))
		result("\n deltaDefocus, defocusMajor, defocusMinor = "+deltaDefocus+" nm, "+defocusMajor+" nm, "+defocusMinor+" nm ")
		//result("\n majorRadius, minorRadius, diff = "+majorRadius+" "+minorRadius+" "+(majorRadius-minorRadius))
		
		//use calcMajor and calcMinor from Fourier transform to represent the length of major and minor radius
		majorRadius1 = calcMajor
		minorRadius1 = calcMinor
		deltaDefocus = (majorRadius1-minorRadius1)*factor
		//displayAt(projectionLines2, 139, 47)  //too many lines to display, SLOW!
		//setWindowSize(projectionLines2, 800, 500)
		
		//displayAt(projectionLines0, 139, 47)  //too many lines to display, SLOW!
		//setWindowSize(projectionLines0, 800, 500)
		
		if (verbose > 1)
		{
			displayAt(ellipseRadius, 139, 47)
			displayAt(ellipseRadiusFFT, 139, 47)
			//displayAt(ellipseRadiusFFTSub, 139, 47)
			displayAt(ellipseRadiusFFTSubModified, 139, 47)
		}
		
		//setWindowSize(ellipseRadius, 800, 400)
		
		Image ret2 := RealImage( "Ellipse_Distortion", 4, 8, 1 )
		ret2[0, 0, 1, 1] = majorDirection //Major axis direction (degree)
		ret2[0, 1, 1, 2] = minorDirection //Minor axis direction (degree)
		ret2[0, 2, 1, 3] = majorRadius1  //Major radius
		ret2[0, 3, 1, 4] = minorRadius1  //Minor radius
		ret2[0, 4, 1, 5] = (majorRadius1-minorRadius1)/(majorRadius1+minorRadius1)/2 //distrotion
		ret2[0, 5, 1, 6] = rangeDeg/polarSamples2 //stepsize
		ret2[0, 6, 1, 7] = mean(ellipseRadius)*factor
		ret2[0, 7, 1, 8] = deltaDefocus  //in nm
		
		result("\n Finish calculatedMinorMajorDirection")

		return ret2


}



Image calculateMinorandMajorRadius(Object self, Number majorDirection, Number minordirection, Number peakPosition, Number s2max, Number shift)
	{
		//oversample = 10
		Number numX, numY
		GetSize(imgS2FFTAmp_FFTAmp, numX, numY)
		result("\n s2max = "+s2max)
		Number rangeDeg = 360
		//find the center of the image
		Number centerX=numX/2, centerY = numY/2
		//determine the half of the smallest dimension
		Number halfMinor = min(numX, numY)/2
		Number polarSamples = rangeDeg * fineFactor3 //use dst30
		Number dr3 = 2*pi() / polarSamples //if rangeDeg = 360, dr = 2*pi()/rangeDeg
		
		//dst30 := createFloatImage( "dst30_Polar_Coordinates", halfMinor, polarSamples )
		//dst30 = warp( imgS2FFTAmp, icol*sin(irow*dr3) + centerX, icol*cos(irow*dr3) + centerY )
		//showimage(dst30)
		//showimage(imgS2FFTAmp)
	
		Number curColMajor = round(majorDirection*(polarSamples/rangeDeg))
		Number curColMinor = round(minorDirection*(polarSamples/rangeDeg))
		result("\n curColMajor = "+curColMajor)
		result("\n curColMinor = "+curColMinor)
		Number majorRadius, minorRadius, temp1, temp2
		
		//Number oversample = 100 //should pad the 1D lines extracted from dst0 (imgS2FFTAmp)
		Number nsamples = halfMinor*(oversample), nx = halfMinor
		//Number peakPositionOversample = peakPosition*(oversample/2)
		Number peakPositionOversample = (peakPosition+shift) * (oversample/2)
		result("\n peakPositionOversample = "+peakPositionOversample)
		result("\n peakPosition*(oversample/2) = "+peakPosition*(oversample/2))
		
		//Number wavelen = 12.2639/sqrt(voltage * 1000.0 + 0.97845 * voltage * voltage)
		
		//Image padded := createFloatImage( "Padded_Projection_One_Line0", nsamples, 1)
		//padded = 0
		//padded[0, 0, 1, floor((oversample/2-0.5)*nx)] = 0
		//padded[0, floor((oversample/2+0.5)*nx), 1, nsamples] = 0
		
		oneLineMinor := dst30[curColMinor, 0, curColMinor+1, halfMinor]
		//ShowImage(oneLineMinor)
		paddedFFTAmpHalf := self.oversample1DArrayShifted(oversample, nx, nsamples, oneLineMinor, shift) 
		//padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = oneLineMinor[icol, irow]
		//Image paddedFFT := RealFFT(padded)
		//Image paddedFFTAmp := (modulus(paddedFFT))
		//showImage(paddedFFTAmp)
		//Image paddedFFTAmpHalf := paddedFFTAmp[0, nsamples/2, 1, nsamples]
		Image mask1D := RealImage( "mask", 4, nsamples/2, 1 )
		mask1D = tert(icol<peakPositionOversample*0.9 || icol>peakPositionOversample*1.1, 0, 1)
		//lowX = peakPositionOversample*0.95
		//highX = peakPositionOversample*1.05
		lowX = peakPositionOversample - 1000
		highX = peakPositionOversample + 1000
		result("\n lowX, highX = "+lowX+" "+highX)
		//mask1D = tert(icol<peakPositionOversample-1000, 0, 1)
		//showImage(mask1D)
		paddedFFTAmpHalf = paddedFFTAmpHalf*mask1D
		max(paddedFFTAmpHalf, minorRadius, temp1)
		Number dfMinor = minorRadius /((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		minorRadius /= (oversample/2) 
		//showImage(paddedFFTAmpHalf)
		
		
		
		//Image padded1 := createFloatImage( "Padded_Projection_One_Line1", nsamples, 1)
		//padded1 = 0
		//padded1[0, 0, 1, floor((oversample/2-0.5)*nx)] = 0
		//padded1[0, floor((oversample/2+0.5)*nx), 1, nsamples] = 0
		oneLineMajor := dst30[curColMajor, 0, curColMajor+1, halfMinor]
		//ShowImage(oneLineMajor)
		paddedFFTAmpHalf1 := self.oversample1DArrayShifted(oversample, nx, nsamples, oneLineMajor, shift) 
		//padded1[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = oneLineMajor[icol, irow]
		//Image paddedFFT1 := RealFFT(padded1)
		//Image paddedFFTAmp1 := (modulus(paddedFFT1))
		//showImage(paddedFFTAmp1)
		//Image paddedFFTAmpHalf1 = paddedFFTAmp1[0, nsamples/2, 1, nsamples] 
		paddedFFTAmpHalf1 = paddedFFTAmpHalf1*mask1D
		max(paddedFFTAmpHalf1, majorRadius, temp1)
		Number dfMajor = majorRadius /((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		majorRadius /= (oversample/2) 
		//showImage(paddedFFTAmpHalf1)
		
		
		result("\n majorRadius after padding = "+majorRadius)
		result("\n minorRadius after padding = "+minorRadius)
		result("\n major defocus = "+dfMajor+ " um")
		result("\n minor defocus = "+dfMinor+" um")
		result("\n delta defocus = "+(dfMajor-dfMinor)+ " um")
		//Number deltaDefocus = (majorRadius - minorRadius)*(oversample/2)/((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		Number deltaDefocus1 = 1*(oversample/2)/((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		result("\n delta defocus of 1 pixel = "+deltaDefocus1+ " um = "+deltaDefocus1*1000+" nm")
		
		image ret3 := RealImage( "Major_and_Minor_Radius", 4, 5, 1 )
		ret3[0, 0, 1, 1] = majorRadius //Major axis 
		ret3[0, 1, 1, 2] = minorRadius //Minor axis
		ret3[0, 2, 1, 3] = (dfMajor - dfMinor) //Minor axis
		ret3[0, 3, 1, 4] = dfMajor
		ret3[0, 4, 1, 5] = dfMinor
		
		return ret3

	}
/*

image calculateMinorandMajorRadius(Object self, Number majorDirection, Number minordirection, Number peakPosition, Number s2max)
	{
		
		Number numX, numY
		GetSize(imgS2FFTAmp_FFTAmp, numX, numY)
		
		Number rangeDeg = 360
		//find the center of the image
		Number centerX=numX/2, centerY = numY/2
		//determine the half of the smallest dimension
		Number halfMinor = min(numX, numY)/2
		Number dr = 2*pi() / rangeDeg //if rangeDeg = 360, dr = 2*pi()/rangeDeg
		Number polarSamples = rangeDeg
	
		Number curColMajor = floor(majorDirection*(polarSamples/rangeDeg))
		Number curColMinor = floor(minorDirection*(polarSamples/rangeDeg))
		Number majorRadius, minorRadius, temp1, temp2
		
		//Number oversample = 100 //should pad the 1D lines extracted from dst0 (imgS2FFTAmp)
		Number nsamples = halfMinor*(oversample), nx = halfMinor
		Number peakPositionOversample = peakPosition*(oversample/2)
		Number wavelen = 12.2639/sqrt(voltage * 1000.0 + 0.97845 * voltage * voltage)
		
		Image padded := createFloatImage( "Padded_Projection_One_Line0", nsamples, 1)
		padded = 0
		padded[0, 0, 1, floor((oversample/2-0.5)*nx)] = 0
		padded[0, floor((oversample/2+0.5)*nx), 1, nsamples] = 0
		
		oneLineMinor := dst0[curColMinor, 0, curColMinor+1, halfMinor]
		
		padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = oneLineMinor[icol, irow]
		Image paddedFFT := RealFFT(padded)
		Image paddedFFTAmp := (modulus(paddedFFT))
		//showImage(paddedFFTAmp)
		Image paddedFFTAmpHalf := paddedFFTAmp[0, nsamples/2, 1, nsamples]
		Image mask1D := RealImage( "mask", 4, nsamples/2, 1 )
		mask1D = tert(icol<peakPositionOversample-1000, 0, 1)
		//showImage(mask1D)
		paddedFFTAmpHalf = paddedFFTAmpHalf*mask1D
		max(paddedFFTAmpHalf, minorRadius, temp1)
		Number dfMinor = minorRadius /((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		minorRadius /= (oversample/2) 
		//showImage(paddedFFTAmpHalf)
		
		
		
		Image padded1 := createFloatImage( "Padded_Projection_One_Line1", nsamples, 1)
		padded1 = 0
		padded1[0, 0, 1, floor((oversample/2-0.5)*nx)] = 0
		padded1[0, floor((oversample/2+0.5)*nx), 1, nsamples] = 0
		oneLineMajor := dst0[curColMajor, 0, curColMajor+1, halfMinor]
		padded1[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = oneLineMajor[icol, irow]
		Image paddedFFT1 := RealFFT(padded1)
		Image paddedFFTAmp1 := (modulus(paddedFFT1))
		//showImage(paddedFFTAmp1)
		Image paddedFFTAmpHalf1 = paddedFFTAmp1[0, nsamples/2, 1, nsamples] 
		paddedFFTAmpHalf1 = paddedFFTAmpHalf1*mask1D
		max(paddedFFTAmpHalf1, majorRadius, temp1)
		Number dfMajor = majorRadius /((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		majorRadius /= (oversample/2) 
		//showImage(paddedFFTAmpHalf1)
		
		
		
		
		result("\n majorRadius after padding = "+majorRadius)
		result("\n minorRadius after padding = "+minorRadius)
		result("\n minor defocus = "+dfMinor+" um")
		result("\n major defocus = "+dfMajor+ " um")
		result("\n delta defocus = "+(dfMajor-dfMinor)+ " um")
		Number deltaDefocus = (majorRadius - minorRadius)*(oversample/2)/((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
		result("\n delta defocus = "+deltaDefocus+ " um")
		
		image ret3 := RealImage( "Major_and_Minor_Radius", 4, 3, 1 )
		ret3[0, 0, 1, 1] = majorRadius //Major axis 
		ret3[0, 1, 1, 2] = minorRadius //Minor axis
		ret3[0, 2, 1, 3] = (dfMajor-dfMinor) //Minor axis
		
		return ret3

	}

*/


void modifiedTrajectory(Object self, Number x, Number y)
{
		Number sizeX, sizeY
		getSize(historyX, sizeX, sizeY)
		//result("\n sizeX = "+sizeX)
		
		for (Number i=sizeX; i>1; i--)
		{
		
			//result("\n i="+i)
			historyX[0, i-1, 1, i] = mean(historyX[0, i-2, 1, i-1])
			historyY[0, i-1, 1, i] = mean(historyY[0, i-2, 1, i-1])
		}
		
		/*
		//transfer the original 3rd element to the 4th element
		historyX[0, 3, 1, 4] = mean(historyX[0, 2, 1, 3])
		historyY[0, 3, 1, 4] = mean(historyY[0, 2, 1, 3])
		
		//transfer the original 2nd element to the 3rd element
		historyX[0, 2, 1, 3] = mean(historyX[0, 1, 1, 2])
		historyY[0, 2, 1, 3] = mean(historyY[0, 1, 1, 2])
		
		//transfer the original 1st element to the 2nd element
		historyX[0, 1, 1, 2] = mean(historyX[0, 0, 1, 1])
		historyY[0, 1, 1, 2] = mean(historyY[0, 0, 1, 1])
		*/
				
		//save the new value (x, y) as the 1st element
		historyX[0, 0, 1, 1] = x
		historyY[0, 0, 1, 1] = y
		
		
		
}


void actOnImage(Object self, Image img0)
	{
	
		String task = ""
		task = "Astigmatism"
		result("\n Your mode is: "+task)
		
		
		



		
		//other important parameters for recording
		//Number beamShiftX, beamShiftY //EMGetBeamShift
		EMGetBeamShift(beamShiftX, beamShiftY)
		//Number calBeamShiftX, calBeamShiftY //EMGetCalibratedBeamShift
		EMGetCalibratedBeamShift(calBeamShiftX, calBeamShiftY)
		//Number beamTiltX, beamTiltY //EMGetBeamTilt
		EMGetBeamTilt(beamTiltX, beamTiltY)
		//Number calBeamTiltX, calBeamTiltY //EMGetCalibratedBeamTilt
		EMGetCalibratedBeamTilt(calBeamTiltX, calBeamTiltY)
		//Number brightness //EMGetBrightness
		brightne = EMGetBrightness(ss)
		//Number getMag //EMGetMagnification
		getMag = EMGetMagnification()
		//Number getFocus, getCalFocus //EMGetFocus, EMGetCalibratedFocus
		getFocus = EMGetFocus()
		getCalFocus = EMGetCalibratedFocus()
		//Number condStigX, condStigY //EMGetCondensorStigmation
		EMGetCondensorStigmation(condStigX, condStigY)
		//Number objStigX, objStigY //EMGetObjectiveStigmation
		EMGetObjectiveStigmation(objStigX, objStigY)
		//Number calObjStigX, calObjStigY //EMGetCalibratedObjectiveStigmation
		EMGetCalibratedObjectiveStigmation(calObjStigX, calObjStigY)
		//Number spotSize //EMGetSpotSize
		spotSize = EMGetSpotSize()
		
		//EMGetObjectiveStigmation(objStigX, objStigY)
		//result("\nobjStigX, objStigY = "+objStigX+", "+objStigY)

		// need to get apixX/apixY from the input image. 
		Number apixX, apixY, numX, numY
		apixX = imageGetDimensionScale( img0, 0 )
		apixY = imageGetDimensionScale( img0, 1 )
		getSize(img0, numX, numY)
		//For DM3, converter=10; for MRC, converter=10000
		apixX = apix
		apixY = apix
		//apixX = apixX*converter
		//apixY = apixY*converter
		//flipVertical(img0)
		
		//perform anisotropic magnification distortion, those parameter will be read from a text file
		//Number aniso = 0.028, theta = 121.3, scale = 1.0  //Mag=14K
		//Number aniso = 0.03, theta = 118.1, scale = 1.0   //Mag=18K
		//Number aniso = 0.026, theta = 118.2, scale = 1.0   //Mag=22.5K
		//Number aniso = 0.026, theta = 114.8, scale = 1.0   //Mag=29K
		//Number aniso = 0.022, theta = 114.9, scale = 1.0   //Mag=37K
		//Number aniso = 0.016, theta = 124.8, scale = 1.0   //Mag=47K
		//Number aniso = 0.004, theta = 150.2, scale = 1.0   //Mag=59K
		
		//Image correctedImg := self.anisotropicScaling(img0, aniso, theta, scale)
		//Image img := correctedImg.imageClone()
		
		Image img := img0.imageClone()
		//ShowImage(img)
		
		
		if (task == "Astigmatism") //input is an image
		{
			
			result("\n actOnImage, your task is "+task)
			Number squareX = min(numX, numY)
			//clip the image to square since the detector is rectangular
			Image imgCloned := img.Slice2(0, 0, 0, 0, squareX, 1, 1, squareX, 1)
			getSize(imgCloned, numX, numY)
			result("\n apix "+apixX+" "+apixY+"\n image size "+numX+"*"+numY) //here numX==numY
			
			
			if (isRGBDataType(img,4))
			  {
			  result("\n The foremost image must not be of type RGB. \n")
			  exit(0)
			  }
			  
			/*if (isComplexDataType(img,8) || isComplexDataType(img,16))
			  {
			  imgS1FFT := imgCloned
			  apixX = 1/(apixX*100)
			  apixY = 1/(apixY*100)
			  }*/
			//else
			  //{
			  imgS1FFT = complexImage("imgS1FFT", 8, numX, numY)
			  convertToFloat(imgCloned)
			  imgS1FFT = RealFFT(imgCloned)
			  deleteImage(imgCloned)
			  
			  imgS1FFTAmp := RealImage( "imgS1_FFT_Amplitude", 4, numX, numY )
			  imgS1FFTAmp = log10(modulus(imgS1FFT))
			  //ShowImage(imgS1FFT) //2D-S1 power spectra
			   
			//  }
			
			
			
			numX = boxsize
			numY = boxsize
			imgS1FFTAmp := RealImage( "imgS1_FFT_Amplitude", 4, numX, numY )
			imgS1FFTAmp = self.gridBoxing(img, boxsize)
			//ShowImage(imgS1FFTAmp)
			
			//prepare to convert 2D-S1 power spectra to 2D-S2 power spectra  
			Number highRes = max(2*apixX, highResolution) //change the default highest resolution limit. Nyquist (2*apix) may include too many background or unhealthy Thon rings, leading to funzzy single ring, especially in high magnification.	
			Number ds, rmax
			ds = 1.0/(apixX*numX)
			rmax = round(1.0/highRes/ds + 0.5)
			if (rmax >= (numX/2)) rmax = numX/2 - 1
			Number s2max = rmax*rmax*ds*ds
			result("\n highRes "+highRes+"\n ds "+ds+"\n rmax "+rmax+"\n s2max "+s2max)
			
			Image ret = self.calculateShiftAndDefocus(imgS1FFTAmp, rmax, s2max, apixX)
			Number recommendMaxShift = mean(ret[0, 0, 1, 1])
			Number meanDefocus = mean(ret[0, 1, 1, 2])
			Number peakPosition = mean(ret[0, 2, 1, 3])
			result("\n Recommended shift = "+recommendMaxShift)
			result("\n Mean defocus = "+meanDefocus+"um")
			result("\n Peak position = "+peakPosition)
			
			
			recommendMaxShift = round(recommendMaxShift*0.75) //75% of the max shift
			//recommendMaxShift = 0
			Image ret2 = self.calculateMinorMajorDirection(imgS2FFTAmp, recommendMaxShift, round(peakPosition), s2max, apixX)
			
			//We only use the major and minor direction returned from function
			Number majorDirection = mean(ret2[0, 0, 1, 1]) //Major axis direction (degree)
			Number minordirection = mean(ret2[0, 1, 1, 2])//Minor axis direction (degree)
			Number majorRadius = mean(ret2[0, 2, 1, 3])
			Number minorRadius = mean(ret2[0, 3, 1, 4])
			//Number distortion = mean(ret2[0, 4, 1, 5]) //Distortion
			Number stepsize = mean(ret2[0, 5, 1, 6]) //stepsize
			Number meanDf = mean(ret2[0, 6, 1, 7]) //stepsize
			Number deltaDefocus = mean(ret2[0, 7, 1, 8]) //astigmatism
			
			result("\n Major Radius after Fourier Shifting: "+majorRadius)
			result("\n Minor Radius after Fourier Shifting: "+minorRadius)
			result("\n stepsize: "+stepsize)
			
			//Image ret3 = self.calculateMinorandMajorRadius(majorDirection, minordirection, round(peakPosition), s2max)
			Image ret3 = self.calculateMinorandMajorRadius(majorDirection, minordirection, round(peakPosition), s2max, recommendMaxShift)
			//majorRadius = mean(ret3[0, 0, 1, 1])
			//minorRadius = mean(ret3[0, 1, 1, 2])
			//Number deltaDefocus = mean(ret3[0, 2, 1, 3])
			//Number dfMajor = mean(ret3[0, 3, 1, 4])
			//Number dfMinor = mean(ret3[0, 4, 1, 5])
			
			
			
			/*
			Number tmpRadius, tmpDirection
			if (majorRadius < minorRadius)
			{
				tmpRadius = majorRadius
				majorRadius = minorRadius
				minorRadius = tmpRadius
				
				tmpDirection = majorDirection
				majorDirection = minordirection
				minordirection = tmpDirection
				
				swap(oneLineMajor, oneLineMinor)
				result("\n Swap major and minor axis!")
			}
			*/
			
			Number distortion = (majorRadius-minorRadius)/((majorRadius+minorRadius)/2)* 100
			//In EM, the default astigAngle is the angle between major axis and +Y. 
			//Here convert it to the angle between major axis and +X. 
			Number astigAngle = (majorDirection+90) < 180 ? (majorDirection+90) : (majorDirection-90)
			
			
			result("\n ***********RESULTS START*******************************************")
			result("\n "+imageName+" apix = "+apix+"  boxsize = "+boxsize+" highRes = "+highResolution+" "+GetTime(1))
			result("\n Mean Defocus (um): "+meanDefocus+" = "+meanDefocus*1000+" nm") //from calculateShiftAndDefocus
			result("\n Mean Defocus (um): "+meanDf+" = "+meanDf+" nm") //from calculateMinorMajorDirection
			//result("\n Astigmatism / Delta defocus (um): "+deltaDefocus+"um = "+deltaDefocus*1000+"nm")
			result("\n Astigmatism / Delta defocus (um): "+deltaDefocus+" nm ")
			result("\n Astigmatism angle to +X (degree): "+astigAngle)
			result("\n Major Axis Direction to +Y (degree): "+majorDirection)
			result("\n Minor Axis Direction to +Y (degree): "+minorDirection)
			result("\n Major Radius in original (pixel): "+majorRadius)
			result("\n Minor Radius in original (pixel): "+minorRadius)
			result("\n Delta Radius in original (pixel): "+(majorRadius - minorRadius))
			result("\n Distortion: "+distortion+"%") //distrotion
			result("\n ***********RESULTS END*********************************************")
			
			wMeanDf = meanDf
			wAstig = deltaDefocus
			wAngle = astigAngle
			
			
			//astigAngle = 180-astigAngle
			
			//plot 2 output figures. 
			Number r_rgb = 0, g_rgb = 1, b_rgb = 0
			Number halfMinor = min(numX, numY)/2
			//Number halfMinorNew = dst3Nx
			//Number offset = halfMinor - halfMinorNew
			img_id = getImageId(ellipseProjectionLines) 
			
			Number nx0, ny0, maxX0, maxY0, maxX1, maxY1
			//Number oversample = 100 //should pad the 1D lines extracted from dst0 (imgS2FFTAmp)
			Number nsamples = halfMinor*(oversample), nx = halfMinor
			nx0 = nsamples/2
			//Number peakPositionOversample = peakPosition*(oversample/2)
			//Image mask1D := RealImage( "mask", 4, nsamples/2, 1 )
			//mask1D = tert(icol<peakPositionOversample-nsamples*0.05, 0, 1)

			
			//showImage(dst20)
			//result(" \n projectionLines0: "+nx0+"*"+ny0)
			if (doesImageExist(img_id)==0 && num_chg==0)
				{
					
					ellipseProjectionLines := createFloatImage("EllipseProjectionLines", nx0, 2)
					ellipseOneLine := createFloatImage("EllipseProjectionAverage", halfMinor, 1)
					//size = 31
					//trajectory := createFloatImage("trajectory", size, size)
					
				}
			if (doesImageExist(img_id)==0 && num_chg>0)
				{
					result("\n Terminating idSTG_beta because ellipseProjectionLines closed. \n")
					exit(0)
				}
				
			if (removePeakLabel) 
				{
					result("\n Remove Peak Label......"+" remove = "+removePeakLabel) 
					lpid.LinePlotImageDisplayRemovePeakLabel( labelPosX )
				}
		
			Number curColMajor = floor(majorDirection)
			Number curColMinor = floor(minorDirection)
			for (Number k=0; k<2; k++)
			{
				r_rgb = (Mod(k,4)==0)?(Abs(r_rgb-1)):r_rgb
				g_rgb = (Mod(k,2)==0)?(Abs(g_rgb-1)):g_rgb
				b_rgb = (Mod(k,1)==0)?(Abs(b_rgb-1)):b_rgb
				if (r_rgb + g_rgb + b_rgb == 3)
					{
					r_rgb = 0
					g_rgb = 0
					b_rgb = 0
					}
				
				
				ellipseOneLine = 0
				
				if (k ==0) 
					{
					
						/*
						Image padded := createFloatImage( "Padded_Projection_One_Line0", nsamples, 1)
						padded = 0
						//padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = ellipseOneLine[icol, irow]
						padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = oneLineMajor[icol, irow]
						Image paddedFFT := RealFFT(padded)
						Image paddedFFTAmp := (modulus(paddedFFT))
						//showImage(paddedFFTAmp)
						Image paddedFFTAmp2 := paddedFFTAmp[0, nsamples/2, 1, nsamples]
						paddedFFTAmp2 = paddedFFTAmp2*mask1D 
						*/
						ellipseProjectionLines[0, 0, 1, nx0] = paddedFFTAmpHalf1
						max(paddedFFTAmpHalf1, maxX0, maxY0)
						//showImage(paddedFFTAmp2)
						
						
					}
				else 
					{
		
						/*
						Image padded := createFloatImage( "Padded_Projection_One_Line0", nsamples, 1)
						padded = 0
						//padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = ellipseOneLine[icol, irow]
						padded[0, (oversample/2-0.5)*nx, 1, (oversample/2+0.5)*nx] = oneLineMinor[icol, irow]
						Image paddedFFT := RealFFT(padded)
						Image paddedFFTAmp := (modulus(paddedFFT))
						//showImage(paddedFFTAmp)
						Image paddedFFTAmp2 := paddedFFTAmp[0, nsamples/2, 1, nsamples] 
						paddedFFTAmp2 = paddedFFTAmp2*mask1D
						*/
						ellipseProjectionLines[1, 0, 2, nx0] = paddedFFTAmpHalf
						max(paddedFFTAmpHalf, maxX1, maxY1)
						
						
					}
					
			
				setDisplayType(ellipseProjectionLines, 4)
				ellipseProjectionLines.ImageSetDimensionUnitString(0, "Defocus (nm)" )
				//ellipseProjectionLines.ImageSetDimensionOrigin(0, -recommendMaxShift*(oversample/2))
				//convert X-axis unit to nm
				Number factor = 1000/((oversample) * s2max * nsamples / (nsamples-1) * 10000 * wavelen)
				ellipseProjectionLines.ImageSetDimensionScale(0, factor)
				ellipseProjectionLines.ImageSetDimensionOrigin(0, -recommendMaxShift*(oversample/2)*factor)
				//ImageSetDimensionScale(ellipseProjectionLines, 0, halfMinor/nx0)
				//ellipseProjectionLines.ImageSetIntensityUnitString( "Intensity" )
				ellipseProjectionLines.ImageSetIntensityUnitString( "Amplitude" )
				
				//linePlotImageDisplay 
				lpid = ellipseProjectionLines.imageGetImageDisplay(0)
				//component lpid = ellipseProjectionLines.imageGetImageDisplay(0)
				lpid.linePlotImageDisplaySetDoAutoSurvey(0,0)
				lpid.linePlotImageDisplaySetGridOn(1)
				lpid.linePlotImageDisplaySetSliceDrawingStyle(k, 1)
				
				//Number offset = nsamples*0.04
				//Number offset = nsamples*0.075 //maybe for high mag
				//Number deltaX = maxX0-maxX1, lowX = min(maxX1, maxX0)-offset, highX = max(maxX1, maxX0)+offset
				lpid.lineplotimagedisplaysetdisplayedchannels(lowX, highX)
				Number minY = min(ellipseProjectionLines), maxY = max(ellipseProjectionLines), deltaY = maxY - minY
				lpid.lineplotimagedisplaysetcontrastlimits(minY - 0.1*deltaY, maxY + 0.25*deltaY)
				//lpid.lineplotimagedisplaysetfontsize(15)
				
				//lpid.lineplotimagedisplaysetlegendshown(1)
				object sliceid0=lpid.imagedisplaygetsliceidbyindex(0)
				lpid.imagedisplaysetslicelabelbyid(sliceid0, "Major Direction "+majorDirection+" degree")
				lpid.lineplotimagedisplaysetsliceLineThickness(sliceid0, 2)
				//lpid.LinePlotImageDisplayAddPeakLabel( maxX0+100, (maxY + 0.025*deltaY), ""+baseN(dfMajor*1000, 10)+"nm" )
			
				object sliceid1=lpid.imagedisplaygetsliceidbyindex(1)
				lpid.imagedisplaysetslicelabelbyid(sliceid1, "Minor Direction "+minorDirection+" degree")
				lpid.lineplotimagedisplaysetsliceLineThickness(sliceid1, 2)
				//lpid.LinePlotImageDisplayAddPeakLabel( maxX1-300, (maxY + 0.025*deltaY), ""+baseN(dfMinor*1000, 10)+"nm" )
				
				Number astig = round(deltaDefocus*10)/10
				//astig = abs(astig)
				Number defocus = round(meanDf*10)/10 + s4Defocus
				if (correctedDefocus) defocus = round(correctedDefocus*10)/10
				
				Number angle = round(astigAngle*10)/10
				String label = "Mean defocus = "+defocus+" nm"+"\nAstigmatism = "+astig+" nm"+"\nAngle = "+angle+" degree"
				labelPosX = (maxX0+maxX1-600)/2
				labelPosY = maxY + 0.2*deltaY
				lpid.LinePlotImageDisplayAddPeakLabel( (maxX0+maxX1-600)/2, maxY + 0.2*deltaY, label )
				
				
				}
				
			removePeakLabel = 1
			
			//Number deltaR = majorRadius-minorRadius
			Number deltaR = deltaDefocus //astigmatism in nm
			Number phiRad = astigAngle/180.0*pi() //in radian
			Number xc = deltaR * cos(phiRad), yc = deltaR * sin(phiRad)
			//Number base2=trunc(log2(deltaR))
			result("\n xc, yc = "+xc+" "+yc)
			
			
			image xVals := [1,1]:
			{
			   { xc }
			}

			image yVals := [1,1]:
			{
			   { yc }
			}
			
			//Number size = 500
			Number sampling = totalAstig/(size/2)  //totalAstig (in nm) ~ 250 pixels
			//size /= sampling
			image xValSampled = round( (xVals) / sampling ) + size/2
			image yValSampled = round( (yVals) / sampling ) + size/2
			result("\n xValSampled, yValSampled = "+sum(xValSampled)+" "+sum(yValSampled))
			
			Number currentX = sum(xValSampled)
			Number currentY = sum(yValSampled)
			
			self.modifiedTrajectory(sum(xValSampled), sum(yValSampled))
			Number dist = distance((sum(xValSampled) - size/2), (sum(yValSampled) - size/2))
			
			if (dist < distance2Origin)
			{
					
					distance2Origin = dist
					bestX = sum(xValSampled)
					bestY = sum(yValSampled)
					result("\n distance = "+dist+" bestX = "+bestX+" bestY = "+bestY)
				
			}
			
			
			
			Number xCtffind3, yCtffind3
			if (astigCtffind3)
			{
				xCtffind3 = round((astigCtffind3 * cos(angCtffind3/180.0*pi())) / sampling) + size/2
				yCtffind3 = round((astigCtffind3 * sin(angCtffind3/180.0*pi())) / sampling) + size/2
				result("\n xCtffind3, yCtffind3 = "+xCtffind3+" "+yCtffind3)
			}
			
			if (doesImageExist(img_id)==0 && num_chg==0)
				{	
					//trajectory := createFloatImage("trajectory", size, size)
					trajectory := RealImage( "Scatter (2D)", 4, size, size )
						
				}
			
			//trajectory := RealImage( "Scatter (2D)", 4, size, size )
			trajectory = 0	
			//trajectory = tert(distance((icol-mean(xValSampled)),(irow-mean(yValSampled))) < 5, 1, 0)
			//trajectory[ xValSampled, yValSampled ] = 1
			
			SetInversionMode(trajectory, 1)
			//showImage(trajectory)
			ImageDisplay imageDisp = trajectory.ImageGetImageDisplay( 0 )
			//Number low, high
			//ImageDisplayGetContrastLimits(imageDisp, low, high)
			//ImageDisplaySetContrastLimits(imageDisp, 0.0, 1.0)
			//flipvertical(trajectory)
			
			//result("\n sum = "+sum(trajectory)+" low = "+low+" high = "+high)
			//result("\n xc, yc, size= "+sum(xValSampled)+" "+sum(yValSampled)+" "+size)
			
			component imgdisp2 = ImageGetImageDisplay(trajectory, 0)
				
			
			if (astigCtffind3)
			{
			//add one circle to represnt the result of ctffind3
			/*
				if (yCtffind3 > size/2) 
					{
						yCtffind3 += size/2
					}
				else
					{
						yCtffind3 += size/2
					}
			*/
				yCtffind3 = yCtffind3 - 2*(yCtffind3-size/2)
				component ovalAnnot5 = NewOvalAnnotation(yCtffind3-6, xCtffind3-6, yCtffind3+6, xCtffind3+6)
				ovalAnnot5.componentsetforegroundcolor(1,0,0)
				ComponentSetFillMode(ovalAnnot5, 1)
				imgdisp2.componentaddchildatend(ovalAnnot5)
			}
			
			if (removeComponents) 
			{
					number numOvals=imgdisp2.componentcountchildrenoftype(6)
					Result("\n numOvals = "+numOvals)
					
					for (Number num=0; num<numOvals; num++)
					{
							component annotid=imgdisp2.componentgetnthchildoftype(6,0)
							annotid.componentremovefromparent()
					}
					
			}
			
			
			
			
			
			Number sizeX, sizeY
			getSize(historyX, sizeX, sizeY)
			
			//plot historical points by different size of the circles
			/*
			Number grayScaleLevel
			Number deltaGrayScale = grayScaleHigh - grayScaleLow
			Number r = 6
			Number nGrayScale2 = 0
			for (Number i=0; i<sizeX; i++)
			{
					Number histX = mean(historyX[0, i, 1, i+1])
					if (histX) 
					{
						nGrayScale2++  //how many non-zero points in the array historyX
					}
					
			}
			result("\n nGrayScale2 = "+nGrayScale2)
			
			for (Number i=0; i<sizeX; i++)
			{
					Number histX = mean(historyX[0, i, 1, i+1]) //the most recent historical point
					Number histY = mean(historyY[0, i, 1, i+1])
					histY = histY - 2*(histY-size/2)
					
					grayScaleLevel = grayScaleLow + (deltaGrayScale/nGrayScale2)*i
					
					
					if (histX) 
					{
							
							//Number r = 2*(sizeX-i)
							//result("\n grayScaleLevel = "+grayScaleLevel)
							component ovalAnnot5 = NewOvalAnnotation(histY-r, histX-r, histY+r, histX+r)
							
							//Number t = (1.0/sizeX)*(sizeX-i)
							//result("\n i = "+i+" histX = "+histX+" histY = "+histY+" t = "+t+" r = "+r)
							//ovalAnnot5.componentsetforegroundcolor(t,t,t)
							ovalAnnot5.componentsetforegroundcolor(grayScaleLevel,grayScaleLevel,grayScaleLevel)
							ComponentSetFillMode(ovalAnnot5, 1)
							imgdisp2.componentaddchildatend(ovalAnnot5)
							
					}
						
			}
			
			
			*/
			
			
			//plot historical points by greyscales of the points
			
			nGrayScale++
			
			Number grayScaleLevel
			Number deltaGrayScale = grayScaleHigh - grayScaleLow
			Number r = 6
			
			for (Number j=0; j<nGrayScale; j++)
			{
				grayScaleLevel = grayScaleLow + (deltaGrayScale/nGrayScale)*j
				//result("\n grayScaleLevel = "+grayScaleLevel)
				
				Number histX2 = mean(historyX[0, j, 1, j+1]) //the most recent historical point
				Number histY2 = mean(historyY[0, j, 1, j+1])
				histY2 = histY2 - 2*(histY2-size/2)
				
				//if (j == (nGrayScale-2))
				//{
				
				component ovalAnnot5 = NewOvalAnnotation(histY2-r, histX2-r, histY2+r, histX2+r)
				ovalAnnot5.componentsetforegroundcolor(grayScaleLevel,grayScaleLevel,grayScaleLevel)
				//ovalAnnot5.componentsetforegroundcolor(0.5,0.5,0.5)
				ComponentSetFillMode(ovalAnnot5, 1)
				imgdisp2.componentaddchildatend(ovalAnnot5)
				//}
				
			}
			
			
			
			result("\n distance = "+dist+" bestX = "+bestX+" bestY = "+bestY)
			Number bestY2 = bestY - 2*(bestY - size/2), bestR = 8
			component ovalAnnot6 = NewOvalAnnotation(bestY2-bestR, bestX-bestR, bestY2+bestR, bestX+bestR)
			ovalAnnot6.componentsetforegroundcolor(1,0,0) //red
			//ComponentSetFillMode(ovalAnnot6, 1)         //fill in red
			imgdisp2.componentaddchildatend(ovalAnnot6)
			
			result("\n currentX = "+currentX+" currentY = "+currentY)
			Number currentY2 = currentY - 2*(currentY - size/2)
			component ovalAnnot7 =  NewOvalAnnotation(currentY2-r, currentX-r, currentY2+r, currentX+r)
			ovalAnnot7.componentsetforegroundcolor(0,0,0)  //black
			ComponentSetFillMode(ovalAnnot7, 1)			   //fill in black
			imgdisp2.componentaddchildatend(ovalAnnot7)
			
			
			flipvertical(trajectory)
			removeComponents = 1
			
			
			
			
			/*
			component ovalAnnot0 = NewOvalAnnotation(0, 0, size, size)
			component ovalAnnot1 = NewOvalAnnotation(size*0.125, size*0.125, size*0.875, size*0.875)
			component ovalAnnot2 = NewOvalAnnotation(size*0.25, size*0.25, size*0.75, size*0.75)
			component ovalAnnot3 = NewOvalAnnotation(size*0.375, size*0.375, size*0.625, size*0.625)
			ovalAnnot0.componentsetdrawingmode(1) // sets the background to off ie lines are not outlined
			ovalAnnot0.componentsetforegroundcolor(0,0,1) 
			ovalAnnot1.componentsetforegroundcolor(0,0,1)
			ovalAnnot2.componentsetforegroundcolor(0,0,1)
			ovalAnnot3.componentsetforegroundcolor(0,0,1)
			imgdisp2.componentaddchildatend(ovalAnnot0)
			imgdisp2.componentaddchildatend(ovalAnnot1)
			imgdisp2.componentaddchildatend(ovalAnnot2)
			imgdisp2.componentaddchildatend(ovalAnnot3)
			*/
			component ovalAnnot0 = NewOvalAnnotation(0, 0, size, size)
			component ovalAnnot1 = NewOvalAnnotation(size*0.1, size*0.1, size*0.9, size*0.9)
			component ovalAnnot2 = NewOvalAnnotation(size*0.2, size*0.2, size*0.8, size*0.8)
			component ovalAnnot3 = NewOvalAnnotation(size*0.3, size*0.3, size*0.7, size*0.7)
			component ovalAnnot4 = NewOvalAnnotation(size*0.4, size*0.4, size*0.6, size*0.6)
			ovalAnnot0.componentsetdrawingmode(1) // sets the background to off ie lines are not outlined
			ovalAnnot0.componentsetforegroundcolor(0.7,0.7,0.7) 
			ovalAnnot1.componentsetforegroundcolor(0.7,0.7,0.7)
			ovalAnnot2.componentsetforegroundcolor(0.7,0.7,0.7)
			ovalAnnot3.componentsetforegroundcolor(0.7,0.7,0.7)
			ovalAnnot4.componentsetforegroundcolor(0.7,0.7,0.7)
			imgdisp2.componentaddchildatend(ovalAnnot0)
			imgdisp2.componentaddchildatend(ovalAnnot1)
			imgdisp2.componentaddchildatend(ovalAnnot2)
			imgdisp2.componentaddchildatend(ovalAnnot3)
			imgdisp2.componentaddchildatend(ovalAnnot4)

			
			component lineAnnot0 = NewLineAnnotation(0, size/2, size, size/2)
			component lineAnnot1 = NewLineAnnotation(size/2, 0, size/2, size)
			lineAnnot0.componentsetforegroundcolor(0.7,0.7,0.7)
			lineAnnot1.componentsetforegroundcolor(0.7,0.7,0.7)
			imgdisp2.componentaddchildatend(lineAnnot0)
			imgdisp2.componentaddchildatend(lineAnnot1)
			/*
			component textAnnot0 = NewTextAnnotation(size/2, 0, ""+totalAstig+"nm", 15)
			component textAnnot1 = NewTextAnnotation(size/2, size/8, ""+totalAstig*0.75+"nm", 15)
			component textAnnot2 = NewTextAnnotation(size/2, size/4, ""+totalAstig*0.5+"nm", 15)
			component textAnnot3 = NewTextAnnotation(size/2, size*3/8, ""+totalAstig*0.25+"nm", 15)
			textAnnot0.componentsetforegroundcolor(1,0,0)
			textAnnot1.componentsetforegroundcolor(1,0,0)
			textAnnot2.componentsetforegroundcolor(1,0,0)
			textAnnot3.componentsetforegroundcolor(1,0,0)
			imgdisp2.componentaddchildatend(textAnnot0)
			imgdisp2.componentaddchildatend(textAnnot1)
			imgdisp2.componentaddchildatend(textAnnot2)
			imgdisp2.componentaddchildatend(textAnnot3)
			*/
			
			component textAnnot0 = NewTextAnnotation(size/2, 0, ""+totalAstig+"nm", 15)
			component textAnnot1 = NewTextAnnotation(size/2, size/10, ""+totalAstig*0.8+"nm", 15)
			component textAnnot2 = NewTextAnnotation(size/2, size*2/10, ""+totalAstig*0.6+"nm", 15)
			component textAnnot3 = NewTextAnnotation(size/2, size*3/10, ""+totalAstig*0.4+"nm", 15)
			component textAnnot4 = NewTextAnnotation(size/2, size*4/10, ""+totalAstig*0.2+"nm", 15)
			textAnnot0.componentsetforegroundcolor(1,0,0)
			textAnnot1.componentsetforegroundcolor(1,0,0)
			textAnnot2.componentsetforegroundcolor(1,0,0)
			textAnnot3.componentsetforegroundcolor(1,0,0)
			textAnnot4.componentsetforegroundcolor(1,0,0)
			imgdisp2.componentaddchildatend(textAnnot0)
			imgdisp2.componentaddchildatend(textAnnot1)
			imgdisp2.componentaddchildatend(textAnnot2)
			imgdisp2.componentaddchildatend(textAnnot3)
			imgdisp2.componentaddchildatend(textAnnot4)
			
			
			//component boxAnnot = NewBoxAnnotation(floor(size/2)-5, floor(size/2)-5, floor(size/2)+5, floor(size/2)+5)
			//boxAnnot.componentsetdrawingmode(1) // sets the background to off ie lines are not outlined
			//boxAnnot.componentsetforegroundcolor(1,0,0) 
			//imgdisp2.componentaddchildatend(boxAnnot)
			ImageDisplaySetContrastLimits(imgdisp2, 0.0, 1.0)
			
			if (verbose > 1) //add an red arrow to represent the direction of astigmatism
			{
				component imgdisp3 = ImageGetImageDisplay(imgS2FFTAmp_FFTAmp, 0)
				Number xp = boxsize/2 + boxsize/2 * cos(phiRad), yp = boxsize/2 - boxsize/2 * sin(phiRad)
				component arrowAnnot3 = NewArrowAnnotation(boxsize/2, boxsize/2, yp, xp)
				arrowAnnot3.componentsetforegroundcolor(1,0,0)
				imgdisp3.componentaddchildatend(arrowAnnot3)
				
				component imgdisp4 = ImageGetImageDisplay(imgS2FFTShifted_FFTAmp, 0)
				component arrowAnnot4 = NewArrowAnnotation(boxsize/2, boxsize/2, yp, xp)
				arrowAnnot4.componentsetforegroundcolor(1,0,0)
				imgdisp4.componentaddchildatend(arrowAnnot4)
				
				component imgdisp5 = ImageGetImageDisplay(imgS1FFTAmpFiltered, 0)
				component arrowAnnot5 = NewArrowAnnotation(boxsize/2, boxsize/2, yp, xp)
				arrowAnnot5.componentsetforegroundcolor(1,0,0)
				imgdisp5.componentaddchildatend(arrowAnnot5)
				
				//showimage(imgS2FFTAmp)
				component imgdisp6 = ImageGetImageDisplay(imgS2FFTAmp0, 0)
				component arrowAnnot6 = NewArrowAnnotation(boxsize/2, boxsize/2, yp, xp)
				arrowAnnot6.componentsetforegroundcolor(1,0,0)
				imgdisp6.componentaddchildatend(arrowAnnot6)
				
				if (clip)
				{
					component imgdisp7 = ImageGetImageDisplay(imgS2FFTAmp_FFTAmp_Clip, 0)
					xp = clip/2 + clip/2 * cos(phiRad)
					yp = clip/2 - clip/2 * sin(phiRad)
					component arrowAnnot7 = NewArrowAnnotation(clip/2, clip/2, yp, xp)
					arrowAnnot7.componentsetforegroundcolor(1,0,0)
					imgdisp7.componentaddchildatend(arrowAnnot7)
					
					component imgdisp8 = ImageGetImageDisplay(imgS1FFTAmpFiltered_Clip, 0)
					component arrowAnnot8 = NewArrowAnnotation(clip/2, clip/2, yp, xp)
					arrowAnnot8.componentsetforegroundcolor(1,0,0)
					imgdisp8.componentaddchildatend(arrowAnnot8)
					
					component imgdisp9 = ImageGetImageDisplay(imgS2FFTAmp0_Clip, 0)
					component arrowAnnot9 = NewArrowAnnotation(clip/2, clip/2, yp, xp)
					arrowAnnot9.componentsetforegroundcolor(1,0,0)
					imgdisp9.componentaddchildatend(arrowAnnot9)
					
				}
				

			}
			
			
			Number imageID
			imageID = getImageId(img)
			
			result("\n imageID = "+imageID+" apix = "+apix+" boxsize = "+boxsize+" highRes = "+highResolution)
			result("\n df = "+wMeanDf+" dfdiff = "+wAstig+" dfang = "+wAngle)
			result("\n nGrayScale = "+nGrayScale)
			result("\n $$$$$$$ No. of points = "+mod(nGrayScale, 5))
			
			
			
			if (saveResults) 
			{
				//WriteFile(fileID, "\n"+imageID+" "+apix+" "+boxsize+" "+clip+" "+highResolution+" "+wMeanDf+" "+wAstig+" "+wAngle)
				WriteFile(fileID, "\n"+imageID+" "+apix+" "+boxsize+" "+clip+" "+highResolution+" "+wMeanDf+" "+wAstig+" "+wAngle+" "+ \
							objStigX+" "+objStigY+" "+calObjStigX+" "+calObjStigY+" "+ \
							getFocus+" "+getCalFocus+" "+condStigX+" "+condStigY+" "+ \
							getMag+" "+brightness+" "+spotSize+" "+ \
							beamShiftX+" "+beamShiftY+" "+calBeamShiftX+" "+calBeamShiftY+" "+ \
							beamTiltX+" "+beamTiltY+" "+calBeamTiltX+" "+calBeamTiltY)
							
							
			}
			
			
			//displayAt(ellipseProjectionLines, 139, 47)
			//setWindowSize(ellipseProjectionLines, 600, 500)
			//SetName(ellipseProjectionLines, imageName)
		
			if (doesImageExist(img_id)==0 && verbose>0)
			 {
				
				
				displayAt(ellipseProjectionLines, 139, 47)
				setWindowSize(ellipseProjectionLines, 500, 400)
				
				displayAt(trajectory, 600, 47)
				//setWindowSize(trajectory, 500, 500)
				
			
				
			  }
			
			updateImage(ellipseProjectionLines)
			//updateImage(trajectory)
			
		
		
		
	
		}
	
		else if (task == "Magnification Distortion")
		{
		
			result("\n actOnImage, you task is "+task)
		}
		
		else
		{
			result("\n No mode is selected. Please select a mode.\n")
			exit(0)
		}
	
	}
	
	
	
	
void actOnChange(object self, number evnt_flg, image img)
   {
		num_chg++
		if (num_chg > 30000)
		  {
		  result("\n Auto exiting due to time-out: re-launch if want to continue. \n")
		  exit(0)
		  }
		  
		Number currentID = getImageId(img)
		//result("\n Current image ID = "+currentID+" Image ID = "+mainImageID)
		/*
		if (currentID != mainImageID) 
			{
				result("\n Current image ID = "+currentID+" Image ID = "+mainImageID+"\n RUN>>>>>>")
				mainImageID = currentID //currentID becomes the new mainImageID
				self.actOnImage(img)
				result("\n I am updating...mean = "+mean(img)+" variance = "+variance(img))
				num_chg = 0
			}
		
		*/
		
		Number repeat
		if (GMS>1) repeat = 5
		else repeat = 8
		
		if (mod(num_chg, repeat)==0)
		  {
			  string evnt_dscr
			  ImageGetEventMap().DeconstructeventFlags( evnt_flg, evnt_dscr )
			  Result("\n "+GetTime(1)+": Image message : " + evnt_dscr + " 0x" + Binary(evnt_flg) + "\n" )
			  self.actOnImage(img)
			  result("\n counter = "+num_chg+" I am updating...mean = "+mean(img)+" variance = "+variance(img))
			  num_chg = 0
		  }
		  
		  
		  
		}
 

} //end of the class idcListenFFT
 


void main()
 {
		object idc_stn
		image img
		string msg_map, task
		number lst_id, img_id
		
		//OKDialog("Welcome to Asitgmatism_MagDistortion.s!"+ \
		//"\nThis script was developed by Rui Yan (yan49@purdue.edu) in Dr. Wen Jiang's group (jiang12@purdue.edu, http://jiang.bio.purdue.edu/) at Purdue University."+ \
		//"\nIt is free for academic study. \nFeel free to contact us if you have any problem.")
		
		//One "trick" used to determine the GMS version since there does not exist a command to get the DM version.
		if (doesFunctionExist("Notes")) 
		{
			result("\n GMS 2.0 or above")
			GMS = 2
			
		}
		else 
		{
			result("\n GMS version is less than 2.0")
			GMS = 1
		}
		  
		if (!img.getFrontImage())
		{
			result("\n No front image found => exiting. \n")
			exit(0)
		}
		
		//open a file to write results
		
		//Number fileID
		if (saveResults)
		{
			
			String filename = "results.txt"
			If (!SaveAsDialog("Save text file as", GetApplicationDirectory(2,0) + "results.txt", filename)) Exit(0)
			Result("\n Selected file path:"+filename)
			fileID = CreateFileForWriting(filename)
			//WriteFile(fileID, "imageID apix boxsize clip highRes df dfdiff dfang")
			WriteFile(fileID, "imageID apix boxsize clip highRes defocus dfdiff dfang objStigX objStigY calObjStigX calObjStigY getFocus getCalFocus condStigX condStigY getMag brightness spotSize beamShiftX beamShiftY calBeamShiftX calBeamShiftY beamTiltX beamTiltY calBeamTiltX calBeamTiltY")
		}
		
		
		imageName = imageGetName(img)	
		img_id = getImageId(img)
		msg_map = "data_changed,data_value_changed:actOnChange"
		idc_stn = alloc(idcListenFFT)
		idc_stn.actOnImage(img)
		lst_id = img.imageAddEventListener(idc_stn, msg_map)
		
		//result("\n imageName = "+imageName+"\n apix = "+apix+"\n boxsize = "+boxsize+"\n highRes = "+highResolution)
		//result("\n df = "+wMeanDf+"\n dfdiff = "+wAstig+"\n dfang = "+wAngle)
		/*
		if (saveResults) 
			{
				WriteFile(fileID, "\n"+imageName+" "+apix+" "+boxsize+" "+clip+" "+highResolution+" "+wMeanDf+" "+wAstig+" "+wAngle)
			}
		*/	
		//exit(0) //for single image test
		while(!shiftDown()) 1==2 //press SHIFT to destroy
		if (saveResults) 
			{
				CloseFile(fileID)
			}
		img.imageRemoveEventListener(lst_id)
		//showImage(historyX)
		//showImage(historyY)
		exit(0)
		
		
		
 }
 
 
 
 

 /*
 void main()
 {
		Object idc_stn
		Image img
		String msg_map, task
		Number lst_id, img_id
		
		String folder, filename
		TagGroup FileList
		number fFiles   = 1
		number fFolders = 2
		
		//One "trick" used to determine the GMS version since there does not exist a command to get the DM version.
		if (doesFunctionExist("Notes")) result("\n GMS 2.0 or above")
		else result("\n GMS version is less than 2.0")
		
		If ( !GetDirectoryDialog( "Select base folder", "", folder ) ) 
			Exit(0)
			
		FileList = GetFilesInDirectory( folder, fFiles + fFolders )
		Result("\n Folder: "+folder)
		
		//open a file to write results
		Number fileID
		if (saveResults)
		{
			
			String filename = "results.txt"
			If (!SaveAsDialog("Save text file as", GetApplicationDirectory(2,0) + "results.txt", filename)) Exit(0)
			Result("\n Selected file path:"+filename)
			fileID = CreateFileForWriting(filename)
			WriteFile(fileID, "imageName apix boxsize clip highRes df dfdiff dfang")
		}
		
		
		number nTags = FileList.TagGroupCountTags()
		Number p1
		for ( number i = 0; i < nTags; i++ )
		{
			TagGroup entryTG
			FileList.TagGroupGetIndexedTagAsTagGroup( i, entryTG )
			if ( entryTG.TagGroupIsValid() )
			{
				string filestr
				if ( entryTG.TagGroupGetTagAsString( "Name", filestr ) )
				{
					
					p1 = find(filestr, ".mrc") //only read MRC files
					//result("\n "+p1)
					if (p1>0) 
						{
							Result( "\n File:" + filestr )
							filename = folder + filestr
							img = OpenImage(filename)
							result("\n Mean = "+mean(img))
							
							imageName = filestr
							//imageName = imageGetName(img)
							//result("\n "+imageName)
							img_id = getImageId(img)
							msg_map = "data_changed,data_value_changed:actOnChange"
							idc_stn = alloc(idcListenFFT)
							idc_stn.actOnImage(img)
							lst_id = img.imageAddEventListener(idc_stn, msg_map)
							//result("\n imageName = "+imageName+"\n apix = "+apix+"\n boxsize = "+boxsize+"\n highRes = "+highResolution)
							//result("\n df = "+wMeanDf+"\n dfdiff = "+wAstig+"\n dfang = "+wAngle)
							if (saveResults) 
							{
								WriteFile(fileID, "\n"+imageName+" "+apix+" "+boxsize+" "+clip+" "+highResolution+" "+wMeanDf+" "+wAstig+" "+wAngle)
							}
							
							deleteImage(img)
							
						}
				}
			}
		}
		
		
		
		
		
		if (saveResults) 
			{
				CloseFile(fileID)
			}
	
 }
 

*/
 

/*

//use for batch processing

void main()
 {
		Object idc_stn
		Image img, img0
		String msg_map, task
		Number lst_id, img_id
		
		//One "trick" used to determine the GMS version since there does not exist a command to get the DM version.
		if (doesFunctionExist("Notes")) result("\n GMS 2.0 or above")
		else result("\n GMS version is less than 2.0")
		
		img.GetFrontImage()
		
		if (!img.getFrontImage())
		{
			result("\n No front image found => exiting. \n")
			exit(0)
		}
		
		Number fileID
		if (saveResults)
		{
			
			String filename = "results.txt"
			If (!SaveAsDialog("Save text file as", GetApplicationDirectory(2,0) + "results.txt", filename)) Exit(0)
			Result("\n Selected file path:"+filename)
			fileID = CreateFileForWriting(filename)
			WriteFile(fileID, "imageName apix boxsize clip highRes df dfdiff dfang")
		}
		
		
		While(img.ImageIsValid())
		{
			imageName = imageGetName(img)
			//result("\n "+imageName)
			img_id = getImageId(img)
			msg_map = "data_changed,data_value_changed:actOnChange"
			idc_stn = alloc(idcListenFFT)
			idc_stn.actOnImage(img)
			lst_id = img.imageAddEventListener(idc_stn, msg_map)
			//result("\n imageName = "+imageName+"\n apix = "+apix+"\n boxsize = "+boxsize+"\n highRes = "+highResolution)
			//result("\n df = "+wMeanDf+"\n dfdiff = "+wAstig+"\n dfang = "+wAngle)
			if (saveResults) 
			{
				WriteFile(fileID, "\n"+imageName+" "+apix+" "+boxsize+" "+clip+" "+highResolution+" "+wMeanDf+" "+wAstig+" "+wAngle)
			}
			
			
			img := FindNextImage(img)
		}
		
		if (saveResults) 
			{
				CloseFile(fileID)
			}
	
 }
 
 */
 
main()