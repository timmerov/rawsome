/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
dump the exposed LibRaw data structures.
**/

#include "dump.h"
#include "log.h"

#include <sstream>

namespace {

std::string quote_or_null(
    const char *str
) {
    if (str == nullptr) {
        return "(nullptr)";
    }
    std::string s("\"");
    s += str;
    s += "\"";
    return s;
}

void dump(
    const std::string& prefix,
    const libraw_image_sizes_t &sizes
) {
    /*
    typedef struct
    {
        ushort raw_height, raw_width, height, width, top_margin, left_margin;
        ushort iheight, iwidth;
        unsigned raw_pitch;
        double pixel_aspect;
        int flip;
        int mask[8][4];
        libraw_raw_crop_t raw_crop;
    } libraw_image_sizes_t;
    */
    LOG(prefix<<".raw_height="<<sizes.raw_height);
    LOG(prefix<<".raw_width="<<sizes.raw_width);
    LOG(prefix<<".height="<<sizes.height);
    LOG(prefix<<".width="<<sizes.width);
    LOG(prefix<<".top_margin="<<sizes.top_margin);
    LOG(prefix<<".left_margin="<<sizes.left_margin);
    LOG(prefix<<".iheight="<<sizes.iheight);
    LOG(prefix<<".iwidth="<<sizes.iwidth);
    LOG(prefix<<".raw_pitch="<<sizes.raw_pitch);
    LOG(prefix<<".pixel_aspect="<<sizes.pixel_aspect);
    LOG(prefix<<".flip="<<sizes.flip);
    std::stringstream ss;
    for (int i = 0; i < 8; ++i) {
        for (int k = 0; k < 4; ++k) {
            ss<<sizes.mask[i][k]<<" ";
        }
    }
    LOG(prefix<<".mask=[ "<<ss.str()<<"]");
    LOG(prefix<<".raw_crop.cleft="<<sizes.raw_crop.cleft);
    LOG(prefix<<".raw_crop.ctop="<<sizes.raw_crop.ctop);
    LOG(prefix<<".raw_crop.cwidth="<<sizes.raw_crop.cwidth);
    LOG(prefix<<".raw_crop.cheight="<<sizes.raw_crop.cheight);
}

void dump(
    const std::string& prefix,
    const libraw_iparams_t &idata
) {
    /*
    typedef struct
    {
        char guard[4];
        char make[64];
        char model[64];
        char software[64];
        unsigned raw_count;
        unsigned dng_version;
        unsigned is_foveon;
        int colors;
        unsigned filters;
        char xtrans[6][6];
        char xtrans_abs[6][6];
        char cdesc[5];
        unsigned xmplen;
        char *xmpdata;
    } libraw_iparams_t;
    */
    std::stringstream ss;
    for (int i = 0; i < 4; ++i) {
        ss<<(int)idata.guard[i]<<" ";
    }
    LOG(prefix<<".guard=[ "<<ss.str()<<"]");
    LOG(prefix<<".make="<<quote_or_null(idata.make));
    LOG(prefix<<".model="<<quote_or_null(idata.model));
    LOG(prefix<<".software="<<quote_or_null(idata.software));
    LOG(prefix<<".raw_count="<<idata.raw_count);
    LOG(prefix<<".dng_version="<<idata.dng_version);
    LOG(prefix<<".is_foveon="<<idata.is_foveon);
    LOG(prefix<<".colors="<<idata.colors);
    LOG(prefix<<".filters="<<idata.filters);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 6; ++k) {
            ss<<(int)idata.xtrans[i][k]<<" ";
        }
    }
    LOG(prefix<<".xtrans=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 6; ++k) {
            ss<<(int)idata.xtrans_abs[i][k]<<" ";
        }
    }
    LOG(prefix<<".xtrans_abs=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<idata.cdesc[i];
    }
    LOG(prefix<<".cdesc='"<<ss.str()<<"'");
    LOG(prefix<<".xmplen="<<idata.xmplen);
    LOG(prefix<<".xmpdata="<<(void*)idata.xmpdata);
}

void dump(
    const std::string& prefix,
    const libraw_lensinfo_t &lens
) {
    /*
    typedef struct
    {
        float MinFocal, MaxFocal, MaxAp4MinFocal, MaxAp4MaxFocal, EXIF_MaxAp;
        char LensMake[128], Lens[128], LensSerial[128], InternalLensSerial[128];
        ushort FocalLengthIn35mmFormat;
        libraw_nikonlens_t nikon;
        libraw_dnglens_t dng;
        libraw_makernotes_lens_t makernotes;
    } libraw_lensinfo_t;
    */
    LOG(prefix<<".MinFocal="<<lens.MinFocal);
    LOG(prefix<<".MaxFocal="<<lens.MaxFocal);
    LOG(prefix<<".MaxAp4MinFocal="<<lens.MaxAp4MinFocal);
    LOG(prefix<<".MaxAp4MaxFocal="<<lens.MaxAp4MaxFocal);
    LOG(prefix<<".EXIF_MaxAp="<<lens.EXIF_MaxAp);
    LOG(prefix<<".LensMake="<<quote_or_null(lens.LensMake));
    LOG(prefix<<".Lens="<<quote_or_null(lens.Lens));
    LOG(prefix<<".LensSerial="<<quote_or_null(lens.LensSerial));
    LOG(prefix<<".InternalLensSerial="<<quote_or_null(lens.InternalLensSerial));
    LOG(prefix<<".FocalLengthIn35mmFormat="<<lens.FocalLengthIn35mmFormat);
    LOG(prefix<<".nikon.NikonEffectiveMaxAp="<<lens.nikon.NikonEffectiveMaxAp);
    LOG(prefix<<".nikon.NikonLensIDNumber="<<(int)lens.nikon.NikonLensIDNumber);
    LOG(prefix<<".nikon.NikonLensFStops="<<(int)lens.nikon.NikonLensFStops);
    LOG(prefix<<".nikon.NikonMCUVersion="<<(int)lens.nikon.NikonMCUVersion);
    LOG(prefix<<".nikon.NikonLensType="<<(int)lens.nikon.NikonLensType);
    LOG(prefix<<".dng.MinFocal="<<(int)lens.dng.MinFocal);
    LOG(prefix<<".dng.MaxFocal="<<(int)lens.dng.MaxFocal);
    LOG(prefix<<".dng.MaxAp4MinFocal="<<(int)lens.dng.MaxAp4MinFocal);
    LOG(prefix<<".dng.MaxAp4MaxFocal="<<(int)lens.dng.MaxAp4MaxFocal);
    /*
    LOG(prefix<<".makernotes.LensID="<<lens.makernotes.LensID);
    LOG(prefix<<".makernotes.Lens="<<quote_or_null(lens.makernotes.Lens));
    LOG(prefix<<".makernotes.LensFormat="<<lens.makernotes.LensFormat);
    LOG(prefix<<".makernotes.LensMount="<<lens.makernotes.LensMount);
    LOG(prefix<<".makernotes.CamID="<<lens.makernotes.CamID);
    LOG(prefix<<".makernotes.CameraFormat="<<lens.makernotes.CameraFormat);
    LOG(prefix<<".makernotes.CameraMount="<<lens.makernotes.CameraMount);
    LOG(prefix<<".makernotes.body="<<quote_or_null(lens.makernotes.body));
    LOG(prefix<<".makernotes.FocalType="<<lens.makernotes.FocalType);
    LOG(prefix<<".makernotes.LensFeatures_pre="<<quote_or_null(lens.makernotes.LensFeatures_pre));
    LOG(prefix<<".makernotes.LensFeatures_suf="<<quote_or_null(lens.makernotes.LensFeatures_suf));
    LOG(prefix<<".makernotes.MinFocal="<<lens.makernotes.MinFocal);
    LOG(prefix<<".makernotes.MaxFocal="<<lens.makernotes.MaxFocal);
    LOG(prefix<<".makernotes.MaxAp4MinFocal="<<lens.makernotes.MaxAp4MinFocal);
    LOG(prefix<<".makernotes.MaxAp4MaxFocal="<<lens.makernotes.MaxAp4MaxFocal);
    LOG(prefix<<".makernotes.MinAp4MinFocal="<<lens.makernotes.MinAp4MinFocal);
    LOG(prefix<<".makernotes.MinAp4MaxFocal="<<lens.makernotes.MinAp4MaxFocal);
    LOG(prefix<<".makernotes.MaxAp="<<lens.makernotes.MaxAp);
    LOG(prefix<<".makernotes.MinAp="<<lens.makernotes.MinAp);
    LOG(prefix<<".makernotes.CurFocal="<<lens.makernotes.CurFocal);
    LOG(prefix<<".makernotes.CurAp="<<lens.makernotes.CurAp);
    LOG(prefix<<".makernotes.MaxAp4CurFocal="<<lens.makernotes.MaxAp4CurFocal);
    LOG(prefix<<".makernotes.MinAp4CurFocal="<<lens.makernotes.MinAp4CurFocal);
    LOG(prefix<<".makernotes.MinFocusDistance="<<lens.makernotes.MinFocusDistance);
    LOG(prefix<<".makernotes.FocusRangeIndex="<<lens.makernotes.FocusRangeIndex);
    LOG(prefix<<".makernotes.LensFStops="<<lens.makernotes.LensFStops);
    LOG(prefix<<".makernotes.TeleconverterID="<<lens.makernotes.TeleconverterID);
    LOG(prefix<<".makernotes.Teleconverter="<<quote_or_null(lens.makernotes.Teleconverter));
    LOG(prefix<<".makernotes.AdapterID="<<lens.makernotes.AdapterID);
    LOG(prefix<<".makernotes.Adapter="<<quote_or_null(lens.makernotes.Adapter));
    LOG(prefix<<".makernotes.AttachmentID="<<lens.makernotes.AttachmentID);
    LOG(prefix<<".makernotes.Attachment="<<quote_or_null(lens.makernotes.Attachment));
    LOG(prefix<<".makernotes.CanonFocalUnits="<<lens.makernotes.CanonFocalUnits);
    LOG(prefix<<".makernotes.FocalLengthIn35mmFormat="<<lens.makernotes.FocalLengthIn35mmFormat);
    */
    LOG(prefix<<".makernotes=not dumped");
}

void dump(
    const std::string& prefix,
    const libraw_shootinginfo_t &shootinginfo
) {
    /*
    typedef struct
    {
        short DriveMode;
        short FocusMode;
        short MeteringMode;
        short AFPoint;
        short ExposureMode;
        short ImageStabilization;
        char BodySerial[64];
        char InternalBodySerial[64]; // this may be PCB or sensor serial, depends on make/model
    } libraw_shootinginfo_t;
    */
    LOG(prefix<<".DriveMode="<<shootinginfo.DriveMode);
    LOG(prefix<<".FocusMode="<<shootinginfo.FocusMode);
    LOG(prefix<<".MeteringMode="<<shootinginfo.MeteringMode);
    LOG(prefix<<".AFPoint="<<shootinginfo.AFPoint);
    LOG(prefix<<".ExposureMode="<<shootinginfo.ExposureMode);
    LOG(prefix<<".ImageStabilization="<<shootinginfo.ImageStabilization);
    LOG(prefix<<".BodySerial="<<quote_or_null(shootinginfo.BodySerial));
    LOG(prefix<<".InternalBodySerial="<<quote_or_null(shootinginfo.InternalBodySerial));
}

void dump(
    const std::string& prefix,
    const libraw_output_params_t &params
) {
    /*
    typedef struct
    {
        unsigned greybox[4];   // -A  x1 y1 x2 y2
        unsigned cropbox[4];   // -B x1 y1 x2 y2
        double aber[4];        // -C
        double gamm[6];        // -g
        float user_mul[4];     // -r mul0 mul1 mul2 mul3
        unsigned shot_select;  // -s
        float bright;          // -b
        float threshold;       //  -n
        int half_size;         // -h
        int four_color_rgb;    // -f
        int highlight;         // -H
        int use_auto_wb;       // -a
        int use_camera_wb;     // -w
        int use_camera_matrix; // +M/-M
        int output_color;      // -o
        char *output_profile;  // -o
        char *camera_profile;  // -p
        char *bad_pixels;      // -P
        char *dark_frame;      // -K
        int output_bps;        // -4
        int output_tiff;       // -T
        int user_flip;         // -t
        int user_qual;         // -q
        int user_black;        // -k
        int user_cblack[4];
        int user_sat; // -S

        int med_passes; // -m
        float auto_bright_thr;
        float adjust_maximum_thr;
        int no_auto_bright;  // -W
        int use_fuji_rotate; // -j
        int green_matching;
        // DCB parameters
        int dcb_iterations;
        int dcb_enhance_fl;
        int fbdd_noiserd;
        int exp_correc;
        float exp_shift;
        float exp_preser;
        // Raw speed
        int use_rawspeed;
        // DNG SDK
        int use_dngsdk;
        // Disable Auto-scale
        int no_auto_scale;
        // Disable intepolation
        int no_interpolation;
        //  int x3f_flags;
        // Sony ARW2 digging mode
        // int sony_arw2_options;
        unsigned raw_processing_options;
        int sony_arw2_posterization_thr;
        // Nikon Coolscan
        float coolscan_nef_gamma;
        char p4shot_order[5];
        // Custom camera list
        char **custom_camera_strings;
    } libraw_output_params_t;
    */
    std::stringstream ss;
    for (int i = 0; i < 4; ++i) {
        ss<<params.greybox[i]<<" ";
    }
    LOG(prefix<<".greybox=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<params.cropbox[i]<<" ";
    }
    LOG(prefix<<".cropbox=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<params.aber[i]<<" ";
    }
    LOG(prefix<<".aber=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 6; ++i) {
        ss<<params.gamm[i]<<" ";
    }
    LOG(prefix<<".gamm=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<params.user_mul[i]<<" ";
    }
    LOG(prefix<<".user_mul=[ "<<ss.str()<<"]");
    LOG(prefix<<".shot_select="<<params.shot_select);
    LOG(prefix<<".bright="<<params.bright);
    LOG(prefix<<".threshold="<<params.threshold);
    LOG(prefix<<".half_size="<<params.half_size);
    LOG(prefix<<".four_color_rgb="<<params.four_color_rgb);
    LOG(prefix<<".highlight="<<params.highlight);
    LOG(prefix<<".use_auto_wb="<<params.use_auto_wb);
    LOG(prefix<<".use_camera_wb="<<params.use_camera_wb);
    LOG(prefix<<".use_camera_matrix="<<params.use_camera_matrix);
    LOG(prefix<<".output_color="<<params.output_color);
    LOG(prefix<<".output_profile="<<quote_or_null(params.output_profile));
    LOG(prefix<<".bad_pixels="<<quote_or_null(params.bad_pixels));
    LOG(prefix<<".dark_frame="<<quote_or_null(params.dark_frame));
    LOG(prefix<<".output_bps="<<params.output_bps);
    LOG(prefix<<".output_tiff="<<params.output_tiff);
    LOG(prefix<<".user_flip="<<params.user_flip);
    LOG(prefix<<".user_qual="<<params.user_qual);
    LOG(prefix<<".user_black="<<params.user_black);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<params.user_cblack[i]<<" ";
    }
    LOG(prefix<<".user_cblack=[ "<<ss.str()<<"]");
    LOG(prefix<<".user_sat="<<params.user_sat);
    LOG(prefix<<".med_passes="<<params.med_passes);
    LOG(prefix<<".auto_bright_thr="<<params.auto_bright_thr);
    LOG(prefix<<".adjust_maximum_thr="<<params.adjust_maximum_thr);
    LOG(prefix<<".no_auto_bright="<<params.no_auto_bright);
    LOG(prefix<<".use_fuji_rotate="<<params.use_fuji_rotate);
    LOG(prefix<<".green_matching="<<params.green_matching);
    LOG(prefix<<".dcb_iterations="<<params.dcb_iterations);
    LOG(prefix<<".dcb_enhance_fl="<<params.dcb_enhance_fl);
    LOG(prefix<<".fbdd_noiserd="<<params.fbdd_noiserd);
    LOG(prefix<<".exp_correc="<<params.exp_correc);
    LOG(prefix<<".exp_shift="<<params.exp_shift);
    LOG(prefix<<".exp_preser="<<params.exp_preser);
    LOG(prefix<<".use_rawspeed="<<params.use_rawspeed);
    LOG(prefix<<".use_dngsdk="<<params.use_dngsdk);
    LOG(prefix<<".no_auto_scale="<<params.no_auto_scale);
    LOG(prefix<<".no_interpolation="<<params.no_interpolation);
    LOG(prefix<<".raw_processing_options="<<params.raw_processing_options);
    LOG(prefix<<".sony_arw2_posterization_thr="<<params.sony_arw2_posterization_thr);
    LOG(prefix<<".coolscan_nef_gamma="<<params.coolscan_nef_gamma);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 5; ++i) {
        ss<<(int)params.p4shot_order[i]<<" ";
    }
    LOG(prefix<<".p4shot_order=[ "<<ss.str()<<"]");
    LOG(prefix<<".custom_camera_strings="<<(void*)params.custom_camera_strings);
}

void dump(
    const std::string& prefix,
    const struct ph1_t &phase_one_data
) {
    /*
    struct ph1_t
    {
        int format, key_off, tag_21a;
        int t_black, split_col, black_col, split_row, black_row;
        float tag_210;
    };
    */
    LOG(prefix<<".format="<<phase_one_data.format);
    LOG(prefix<<".key_off="<<phase_one_data.key_off);
    LOG(prefix<<".tag_21a="<<phase_one_data.tag_21a);
    LOG(prefix<<".t_black="<<phase_one_data.t_black);
    LOG(prefix<<".split_col="<<phase_one_data.split_col);
    LOG(prefix<<".black_col="<<phase_one_data.black_col);
    LOG(prefix<<".split_row="<<phase_one_data.split_row);
    LOG(prefix<<".black_row="<<phase_one_data.black_row);
    LOG(prefix<<".tag_210="<<phase_one_data.tag_210);
}

void dump(
    const std::string& prefix,
    const libraw_dng_color_t &dng_color
) {
    LOG(prefix<<".parsedfields="<<dng_color.parsedfields);
    LOG(prefix<<".illuminant="<<dng_color.illuminant);
    std::stringstream ss;
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
            ss<<dng_color.calibration[i][k]<<" ";
        }
    }
    LOG(prefix<<".calibration=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 3; ++k) {
            ss<<dng_color.colormatrix[i][k]<<" ";
        }
    }
    LOG(prefix<<".colormatrix=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 4; ++k) {
            ss<<dng_color.forwardmatrix[i][k]<<" ";
        }
    }
    LOG(prefix<<".forwardmatrix=[ "<<ss.str()<<"]");
}

void dump(
    const std::string& prefix,
    const libraw_dng_levels_t &dng_levels
) {
    /*
    typedef struct
    {
        unsigned parsedfields;
        unsigned dng_cblack[4102];
        unsigned dng_black;
        unsigned dng_whitelevel[4];
        unsigned default_crop[4]; // Origin and size
        unsigned preview_colorspace;
        float analogbalance[4];
    } libraw_dng_levels_t;
    */
    LOG(prefix<<".parsedfields="<<dng_levels.parsedfields);
    LOG(prefix<<".dng_cblack=not dumped");
    LOG(prefix<<".dng_black="<<dng_levels.dng_black);
    std::stringstream ss;
    for (int i = 0; i < 4; ++i) {
        ss<<dng_levels.dng_whitelevel[i]<<" ";
    }
    LOG(prefix<<".dng_whitelevel=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<dng_levels.default_crop[i]<<" ";
    }
    LOG(prefix<<".default_crop=[ "<<ss.str()<<"]");
    LOG(prefix<<".preview_colorspace="<<dng_levels.preview_colorspace);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<dng_levels.analogbalance[i]<<" ";
    }
    LOG(prefix<<".analogbalance=[ "<<ss.str()<<"]");
}

void dump(
    const std::string& prefix,
    const libraw_P1_color_t &P1_color
) {
    /*
    typedef struct
    {
        float romm_cam[9];
    } libraw_P1_color_t;
    */
    std::stringstream ss;
    for (int i = 0; i < 9; ++i) {
        ss<<P1_color.romm_cam[i]<<" ";
    }
    LOG(prefix<<".calibration=[ "<<ss.str()<<"]");
}

void dump(
    const std::string& prefix,
    const libraw_colordata_t &color
) {
    /*
    typedef struct
    {
        ushort curve[0x10000];
        unsigned cblack[4102];
        unsigned black;
        unsigned data_maximum;
        unsigned maximum;
        long linear_max[4];
        float fmaximum;
        float fnorm;
        ushort white[8][8];
        float cam_mul[4];
        float pre_mul[4];
        float cmatrix[3][4];
        float ccm[3][4];
        float rgb_cam[3][4];
        float cam_xyz[4][3];
        struct ph1_t phase_one_data;
        float flash_used;
        float canon_ev;
        char model2[64];
        char UniqueCameraModel[64];
        char LocalizedCameraModel[64];
        void *profile;
        unsigned profile_length;
        unsigned black_stat[8];
        libraw_dng_color_t dng_color[2];
        libraw_dng_levels_t dng_levels;
        float baseline_exposure;
        int WB_Coeffs[256][4];    // R, G1, B, G2 coeffs
        float WBCT_Coeffs[64][5]; // CCT, than R, G1, B, G2 coeffs
        libraw_P1_color_t P1_color[2];
    } libraw_colordata_t;
    */
    std::stringstream ss;
    for (int i = 0; i < 0x10000; i += 0x100) {
        ss<<color.curve[i]<<" ";
    }
    LOG(prefix<<".curve=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 6; ++i) {
        ss<<color.cblack[i]<<" ";
    }
    LOG(prefix<<".cblack=[ "<<ss.str()<<"... ]");
    LOG(prefix<<".black="<<color.black);
    LOG(prefix<<".data_maximum="<<color.data_maximum);
    LOG(prefix<<".maximum="<<color.maximum);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 5; ++i) {
        ss<<color.linear_max[i]<<" ";
    }
    LOG(prefix<<".linear_max=[ "<<ss.str()<<"]");
    LOG(prefix<<".fmaximum="<<color.fmaximum);
    LOG(prefix<<".fnorm="<<color.fnorm);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 8; ++i) {
        for (int k = 0; k < 8; ++k) {
            ss<<color.white[i][k]<<" ";
        }
    }
    LOG(prefix<<".white=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<color.cam_mul[i]<<" ";
    }
    LOG(prefix<<".cam_mul=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<color.pre_mul[i]<<" ";
    }
    LOG(prefix<<".pre_mul=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 4; ++k) {
            ss<<color.cmatrix[i][k]<<" ";
        }
    }
    LOG(prefix<<".cmatrix=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 4; ++k) {
            ss<<color.ccm[i][k]<<" ";
        }
    }
    LOG(prefix<<".ccm=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 4; ++k) {
            ss<<color.rgb_cam[i][k]<<" ";
        }
    }
    LOG(prefix<<".rgb_cam=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 3; ++k) {
            ss<<color.cam_xyz[i][k]<<" ";
        }
    }
    LOG(prefix<<".cam_xyz=[ "<<ss.str()<<"]");
    dump(prefix+".phase_one_data", color.phase_one_data);
    LOG(prefix<<".flash_used="<<color.flash_used);
    LOG(prefix<<".canon_ev="<<color.canon_ev);
    LOG(prefix<<".model2="<<quote_or_null(color.model2));
    LOG(prefix<<".UniqueCameraModel="<<quote_or_null(color.UniqueCameraModel));
    LOG(prefix<<".LocalizedCameraModel="<<quote_or_null(color.LocalizedCameraModel));
    LOG(prefix<<".profile="<<(void*)color.profile);
    LOG(prefix<<".profile_length="<<color.profile_length);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 8; ++i) {
        ss<<color.black_stat[i]<<" ";
    }
    LOG(prefix<<".black_stat=[ "<<ss.str()<<"]");
    dump(prefix+".dng_color[0]", color.dng_color[0]);
    dump(prefix+".dng_color[1]", color.dng_color[1]);
    dump(prefix+".dng_levels", color.dng_levels);
    LOG(prefix<<".baseline_exposure="<<color.baseline_exposure);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 256; ++i) {
        for (int k = 0; k < 4; ++k) {
            ss<<color.WB_Coeffs[i][k]<<" ";
        }
    }
    LOG(prefix<<".WB_Coeffs=[ "<<ss.str()<<"]");
    LOG(prefix<<".WBCT_Coeffs=not dumped");
    dump(prefix+".P1_color[0]", color.P1_color[0]);
    dump(prefix+".P1_color[1]", color.P1_color[1]);
}

void dump(
    const std::string& prefix,
    const libraw_gps_info_t &parsed_gps
) {
    /*
    typedef struct
    {
        float latitude[3];     // Deg,min,sec
        float longtitude[3];   // Deg,min,sec
        float gpstimestamp[3]; // Deg,min,sec
        float altitude;
        char altref, latref, longref, gpsstatus;
        char gpsparsed;
    } libraw_gps_info_t;
    */
    std::stringstream ss;
    for (int i = 0; i < 3; ++i) {
        ss<<parsed_gps.latitude[i]<<" ";
    }
    LOG(prefix<<".latitude=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        ss<<parsed_gps.longtitude[i]<<" ";
    }
    LOG(prefix<<".longtitude=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        ss<<parsed_gps.gpstimestamp[i]<<" ";
    }
    LOG(prefix<<".gpstimestamp=[ "<<ss.str()<<"]");
    LOG(prefix<<".altitude="<<parsed_gps.altitude);
    LOG(prefix<<".altref="<<(int)parsed_gps.altref);
    LOG(prefix<<".latref="<<(int)parsed_gps.latref);
    LOG(prefix<<".longref="<<(int)parsed_gps.longref);
    LOG(prefix<<".gpsstatus="<<(int)parsed_gps.gpsstatus);
    LOG(prefix<<".gpsparsed="<<(int)parsed_gps.gpsparsed);
}

void dump(
    const std::string& prefix,
    const libraw_imgother_t &other
) {
    /*
    typedef struct
    {
        float iso_speed;
        float shutter;
        float aperture;
        float focal_len;
        time_t timestamp;
        unsigned shot_order;
        unsigned gpsdata[32];
        libraw_gps_info_t parsed_gps;
        char desc[512], artist[64];
        float FlashEC;
        float FlashGN;
        float CameraTemperature;
        float SensorTemperature;
        float SensorTemperature2;
        float LensTemperature;
        float AmbientTemperature;
        float BatteryTemperature;
        float exifAmbientTemperature;
        float exifHumidity;
        float exifPressure;
        float exifWaterDepth;
        float exifAcceleration;
        float exifCameraElevationAngle;
        float real_ISO;
    } libraw_imgother_t;
    */
    LOG(prefix<<".iso_speed="<<other.iso_speed);
    LOG(prefix<<".shutter="<<other.shutter);
    LOG(prefix<<".aperture="<<other.aperture);
    LOG(prefix<<".focal_len="<<other.focal_len);
    LOG(prefix<<".timestamp="<<other.timestamp);
    LOG(prefix<<".shot_order="<<other.shot_order);
    std::stringstream ss;
    for (int i = 0; i < 32; ++i) {
        ss<<other.gpsdata[i]<<" ";
    }
    LOG(prefix<<".gpsdata=[ "<<ss.str()<<"]");
    dump(prefix+".parsed_gps", other.parsed_gps);
    LOG(prefix<<".desc="<<quote_or_null(other.desc));
    LOG(prefix<<".artist="<<quote_or_null(other.artist));
    LOG(prefix<<".FlashEC="<<other.FlashEC);
    LOG(prefix<<".FlashGN="<<other.FlashGN);
    LOG(prefix<<".CameraTemperature="<<other.CameraTemperature);
    LOG(prefix<<".SensorTemperature="<<other.SensorTemperature);
    LOG(prefix<<".SensorTemperature2="<<other.SensorTemperature2);
    LOG(prefix<<".LensTemperature="<<other.LensTemperature);
    LOG(prefix<<".AmbientTemperature="<<other.AmbientTemperature);
    LOG(prefix<<".BatteryTemperature="<<other.BatteryTemperature);
    LOG(prefix<<".exifAmbientTemperature="<<other.exifAmbientTemperature);
    LOG(prefix<<".exifHumidity="<<other.exifHumidity);
    LOG(prefix<<".exifPressure="<<other.exifPressure);
    LOG(prefix<<".exifWaterDepth="<<other.exifWaterDepth);
    LOG(prefix<<".exifAcceleration="<<other.exifAcceleration);
    LOG(prefix<<".exifCameraElevationAngle="<<other.exifCameraElevationAngle);
    LOG(prefix<<".real_ISO="<<other.real_ISO);
}

void dump(
    const std::string& prefix,
    const libraw_thumbnail_t &thumbnail
) {
    /*
    typedef struct
    {
        enum LibRaw_thumbnail_formats tformat;
        ushort twidth, theight;
        unsigned tlength;
        int tcolors;
        char *thumb;
    } libraw_thumbnail_t;
    */
    LOG(prefix<<".tformat="<<(int)thumbnail.tformat);
    LOG(prefix<<".twidth="<<thumbnail.twidth);
    LOG(prefix<<".theight="<<thumbnail.theight);
    LOG(prefix<<".tlength="<<thumbnail.tlength);
    LOG(prefix<<".tcolors="<<thumbnail.tcolors);
    LOG(prefix<<".thumb="<<(void*)thumbnail.thumb);
}

void dump(
    const std::string& prefix,
    const libraw_internal_output_params_t &ioparams
) {
    /*
    typedef struct
    {
        unsigned mix_green;
        unsigned raw_color;
        unsigned zero_is_bad;
        ushort shrink;
        ushort fuji_width;
    } libraw_internal_output_params_t;
    */
    LOG(prefix<<".mix_green="<<ioparams.mix_green);
    LOG(prefix<<".raw_color="<<ioparams.raw_color);
    LOG(prefix<<".zero_is_bad="<<ioparams.zero_is_bad);
    LOG(prefix<<".shrink="<<ioparams.shrink);
    LOG(prefix<<".fuji_width="<<ioparams.fuji_width);
}

void dump(
    const std::string& prefix,
    const libraw_rawdata_t &rawdata
) {
    /*
    typedef struct
    {
        // really allocated bitmap
        void *raw_alloc;
        // alias to single_channel variant
        ushort *raw_image;
        // alias to 4-channel variant
        ushort (*color4_image)[4];
        // alias to 3-color variand decoded by RawSpeed
        ushort (*color3_image)[3];
        // float bayer
        float *float_image;
        // float 3-component
        float (*float3_image)[3];
        // float 4-component
        float (*float4_image)[4];

        // Phase One black level data;
        short (*ph1_cblack)[2];
        short (*ph1_rblack)[2];
        // save color and sizes here, too....
        libraw_iparams_t iparams;
        libraw_image_sizes_t sizes;
        libraw_internal_output_params_t ioparams;
        libraw_colordata_t color;
    } libraw_rawdata_t;
    */
    LOG(prefix<<".raw_alloc="<<(void*)rawdata.raw_alloc);
    LOG(prefix<<".raw_image="<<(void*)rawdata.raw_image);
    std::stringstream ss;
    for (int i = 0; i < 4; ++i) {
        ss<<(void*)rawdata.color4_image[i]<<" ";
    }
    LOG(prefix<<".color4_image=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        ss<<(void*)rawdata.color3_image[i]<<" ";
    }
    LOG(prefix<<".color3_image=[ "<<ss.str()<<"]");
    LOG(prefix<<".float_image="<<(void*)rawdata.float_image);
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 3; ++i) {
        ss<<(void*)rawdata.float3_image[i]<<" ";
    }
    LOG(prefix<<".float3_image=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 4; ++i) {
        ss<<(void*)rawdata.float4_image[i]<<" ";
    }
    LOG(prefix<<".float4_image=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 2; ++i) {
        ss<<(void*)rawdata.ph1_cblack[i]<<" ";
    }
    LOG(prefix<<".ph1_cblack=[ "<<ss.str()<<"]");
    ss.clear();
    ss.str(std::string());
    for (int i = 0; i < 2; ++i) {
        ss<<(void*)rawdata.ph1_rblack[i]<<" ";
    }
    LOG(prefix<<".ph1_rblack=[ "<<ss.str()<<"]");
    dump(prefix+".iparams", rawdata.iparams);
    dump(prefix+".sizes", rawdata.sizes);
    dump(prefix+".ioparams", rawdata.ioparams);
    dump(prefix+".color", rawdata.color);
}

void dump(
    const std::string& prefix,
    const libraw_data_t &imgdata
) {
    /*
    typedef struct
    {
        ushort (*image)[4];
        libraw_image_sizes_t sizes;
        libraw_iparams_t idata;
        libraw_lensinfo_t lens;
        libraw_makernotes_t makernotes;
        libraw_shootinginfo_t shootinginfo;
        libraw_output_params_t params;
        unsigned int progress_flags;
        unsigned int process_warnings;
        libraw_colordata_t color;
        libraw_imgother_t other;
        libraw_thumbnail_t thumbnail;
        libraw_rawdata_t rawdata;
        void *parent_class;
    } libraw_data_t;
    */
    LOG(prefix<<".image[0]="<<(void*)imgdata.image[0]);
    LOG(prefix<<".image[1]="<<(void*)imgdata.image[1]);
    LOG(prefix<<".image[2]="<<(void*)imgdata.image[2]);
    LOG(prefix<<".image[3]="<<(void*)imgdata.image[3]);
    dump(prefix+".sizes", imgdata.sizes);
    dump(prefix+".idata", imgdata.idata);
    dump(prefix+".lens", imgdata.lens);
    LOG(prefix<<".makernotes=not dumped");
    dump(prefix+".shootinginfo", imgdata.shootinginfo);
    dump(prefix+".params", imgdata.params);
    LOG(prefix<<".progress_flags="<<imgdata.progress_flags);
    LOG(prefix<<".process_warnings="<<imgdata.process_warnings);
    dump(prefix+".color", imgdata.color);
    dump(prefix+".other", imgdata.other);
    dump(prefix+".thumbnail", imgdata.thumbnail);
    dump(prefix+".rawdata", imgdata.rawdata);
    LOG(prefix<<".parent_class="<<(void*)imgdata.parent_class);
}

}

void dump(
    const LibRaw &raw_image
) {
    LOG("dumping raw_image:");
    dump("imgdata", raw_image.imgdata);
}
