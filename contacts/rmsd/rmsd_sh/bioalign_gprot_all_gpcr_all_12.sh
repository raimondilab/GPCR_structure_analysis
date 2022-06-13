#!/bin/bash
./bio_align_fit_n_rmsd.py -refe ./cifs/7wvy.cif -mobi ./cifs/7xji.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,144,145,146,147,148,149,150,151,152,153,154,155,157,158,159,203,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,257,258,259,260,261,262,263,264,265,266,267,268,269,321,322,323,324,332 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7wvy_7xji_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9a.cif -mobi ./cifs/7x9b.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,213,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278 -m_fit_sele 56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,167,168,169,170,171,172,173,174,175,176,177,178,180,181,182,221,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 > 7x9a_7x9b_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9a.cif -mobi ./cifs/7x9c.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,213,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278 -m_fit_sele 47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,216,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 > 7x9a_7x9c_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9a.cif -mobi ./cifs/7xjh.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,213,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7x9a_7xjh_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9a.cif -mobi ./cifs/7xji.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,213,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7x9a_7xji_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9b.cif -mobi ./cifs/7x9c.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,167,168,169,170,171,172,173,174,175,176,177,178,180,181,182,221,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283 -m_fit_sele 47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,216,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 > 7x9b_7x9c_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9b.cif -mobi ./cifs/7xjh.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,167,168,169,170,171,172,173,174,175,176,177,178,180,181,182,221,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7x9b_7xjh_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9b.cif -mobi ./cifs/7xji.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,167,168,169,170,171,172,173,174,175,176,177,178,180,181,182,221,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7x9b_7xji_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9c.cif -mobi ./cifs/7xjh.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,216,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7x9c_7xjh_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7x9c.cif -mobi ./cifs/7xji.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,216,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,256,257,258,259,260,261,262,263,264,265,266,267,268,320,321,322,323,331 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7x9c_7xji_fit.txt
./bio_align_fit_n_rmsd.py -refe ./cifs/7xjh.cif -mobi ./cifs/7xji.cif -r_fit_chain R -m_fit_chain R -r_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -m_fit_sele 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,156,157,158,159,160,161,162,163,164,165,166,167,169,171,172,206,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307 -r_rms_chain A -m_rms_chain A -r_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 -m_rms_sele 209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,269,270,271,272,273,274,275,276,277,278,279,280,281,350,351,352,353,361 > 7xjh_7xji_fit.txt
