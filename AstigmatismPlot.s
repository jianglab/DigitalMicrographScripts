// $BACKGROUND$
// Plots 1/px vs. mean log10(norm(I)) for a number of angles
// if front image is complex, it is assumed to be the FFT
// otherwise, it's FFT is taken here
// Written by: ilker donmez around Nov 1, 10
// Benefited from various scripts at: http://www.felmi-zfe.tugraz.at/dm_scripts/welcome.html

class idcListenFFT
  {
  number num_deg, num_chg, apix_x, apix_y, num_x, num_y
  number cent_x, cent_y, inc_deg, hlf_mnr, num_wdg, cur_wdg
  number r_rgb, g_rgb, b_rgb, ml2v, img_id, img_type
  image tmp, plr_img, prj_lns, rot_avg, tmp_frn, lns_sub
  compleximage img_frn, img_fft
  // currently not used: apix_x, apix_y, img_type
  
  idcListenFFT(object self)
    {
    num_chg = 0
    result("\n object 0x"+self.scriptObjectGetID().hex()+" created.")
    }
  
  ~idcListenFFT(object self)
    {
    result("\n object 0x"+self.scriptObjectGetID().hex()+" destroyed.")
    }
  
  void actOnImage(object self, image img)
    {
    num_deg = 180
    num_wdg = 4
    r_rgb = 0
    g_rgb = 1
    b_rgb = 0
    tmp_frn = imageClone(img)
    apix_x = imageGetDimensionScale( img, 0)
    apix_y = imageGetDimensionScale( img, 1)
    // need to get apix_x/y from orginal img
    getSize(tmp_frn, num_x, num_y)
    // img_type = imageGetDataType(tmp_frn)
    ml2v = mod(log2(num_x),1)
    if (num_x != num_y || ml2v != 0)
      {
      result("\n Image must be square and with 2^int dimensions. \n")
      exit(0)
      }
    if (isRGBDataType(img,4))
      {
      result("\n The foremost image must not be of type RGB. \n")
      exit(0)
      }
    if (isComplexDataType(img,8) || isComplexDataType(img,16))
      {
      img_frn := tmp_frn
      apix_x = 1/(apix_x*100)
      apix_y = 1/(apix_y*100)
      }
    else
      {
      img_fft = complexImage("img_fft", 8, num_x, num_y)
      convertToFloat(tmp_frn)
      img_fft = realFft(tmp_frn)
      deleteImage(tmp_frn)
      img_frn := img_fft
      apix_x = apix_x*10
      apix_y = apix_y*10
      }
    convertToComplex(img_frn)
    tmp := log10(modulus(img_frn))
    hlf_mnr = min( num_x, num_y )/2
    cent_x = num_x / 2
    cent_y = num_y / 2
    plr_img := createFloatImage( "plr_img", hlf_mnr, num_deg )
    inc_deg = 1*pi() / num_deg
    plr_img = warp( tmp, icol*sin(irow*inc_deg) + cent_x, icol*cos(irow*inc_deg) + cent_y )
    img_id = getImageId(prj_lns)
    if (doesImageExist(img_id)==0 && num_chg==0)
      {
      // deleteImage(prj_lns)
      prj_lns := createFloatImage( "prj_lns", hlf_mnr, num_wdg)
      }
    if (doesImageExist(img_id)==0 && num_chg>0)
      {
      result("\n Terminating idSTG_beta because prj_lns closed. \n")
      // self.~idcListenFFT()
      exit(0)
      }
    prj_lns = 0
    rot_avg := createFloatImage("rot_avg", hlf_mnr, num_deg/num_wdg )
    lns_sub = createFloatImage("lns_sub", hlf_mnr, 1)
    
    for (cur_wdg=0; cur_wdg<num_wdg; cur_wdg++)
      {
      r_rgb = (Mod(cur_wdg,4)==0)?(Abs(r_rgb-1)):r_rgb
      g_rgb = (Mod(cur_wdg,2)==0)?(Abs(g_rgb-1)):g_rgb
      b_rgb = (Mod(cur_wdg,1)==0)?(Abs(b_rgb-1)):b_rgb
      if (r_rgb + g_rgb + b_rgb == 3)
        {
        r_rgb = 0
        g_rgb = 0
        b_rgb = 0
        }
      rot_avg = 0
      rot_avg = plr_img[cur_wdg*num_deg/num_wdg, 0, (cur_wdg+1)*num_deg/num_wdg, hlf_mnr]
      lns_sub = 0
      lns_sub[icol, irow] += rot_avg
      prj_lns[cur_wdg,0,cur_wdg+1,hlf_mnr] = lns_sub
      prj_lns[cur_wdg,0,cur_wdg+1,hlf_mnr] /= num_deg/num_wdg
      // setName(prj_lns, "prj_lns")
      // this is not necessary
      setDisplayType(prj_lns, 4)
      prj_lns.ImageSetDimensionUnitString(0, "1/A-depends on image info in DM" )
      SetScale(prj_lns, 1/(apix_x*(num_x-1)), 1/(apix_y*(num_y-1)))
      // prj_lns.imageSetDimensionUnitString(0, "1/px" )
      // setScale(prj_lns, 1/(num_x-1), 1/(num_x-1))
      linePlotImageDisplay lpid = prj_lns.imageGetImageDisplay(0);
      lpid.linePlotImageDisplaySetDoAutoSurvey(1,1);
      lpid.linePlotImageDisplaySetGridOn(0);
      lpid.linePlotImageDisplaySetSliceDrawingStyle(cur_wdg, 1);
      lpid.linePlotImageDisplaySetSliceComponentColor(cur_wdg, 0, r_rgb, g_rgb, b_rgb);
      }
    if (doesImageExist(img_id)==0)
      {
      displayAt(prj_lns, 139, 47);
      setWindowSize(prj_lns, 800, 325);
      }
    updateImage(prj_lns)
    }
  
  void actOnChange(object self, number evnt_flg, image img)
    {
    num_chg++
    if (num_chg > 30000)
      {
      Result("\n Auto exiting due to time-out: re-launch if want to continue. \n")
      exit(0)
      }
    if (mod(num_chg,1)==0)
      {
      //string evnt_dscr
      //ImageGetEventMap().DeconstructeventFlags( evnt_flg, evnt_dscr )
      //Result(GetTime(1)+": Image message : " + evnt_dscr + " 0x" + Binary(evnt_flg) + "\n" )
      self.actOnImage(img)
      }
    }
  }


void main()
  {
  object idc_stn
  image img
  string msg_map
  number lst_id, img_id
  
  if (!img.getFrontImage())
    {
    result("\n No front image found => exiting. \n")
    exit(0)
    }
  img_id = getImageId(img)
  msg_map = "data_changed,data_value_changed:actOnChange"
  idc_stn = alloc(idcListenFFT)
  idc_stn.actOnImage(img)
  lst_id = img.imageAddEventListener(idc_stn, msg_map)
  exit(0)
  while(!shiftDown()) 1==2
  img.imageRemoveEventListener(lst_id)
  exit(0)
  // last three lines not used
  }


main()
