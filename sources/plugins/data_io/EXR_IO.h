#ifndef EXR_IO_H_
#define EXR_IO_H_
/*
 * Author: Cyril Soler
 */


#include <stdexcept>
#include <cassert>

#include <ImfRgbaFile.h>
#include <ImfArray.h>

// XXX: This header is not installed as of version 2.2.0 of OpenEXR, see
// <https://lists.nongnu.org/archive/html/openexr-devel/2016-06/msg00001.html>
// and <https://github.com/openexr/openexr/pull/184>.
#include <ImfStdIO.h>

template<typename FType>
class t_EXR_IO
{
	public:
		static bool LoadEXR(std::istream& input, int& W,int& H, FType *& pix,int nC=3) 
		{
      // XXX: OpenEXR implements its own IStream and OStream classes, but
      // they are completely independent from those of libstdc++.  The
      // closest thing it has is 'StdIFStream', hence this hack.
      std::ifstream* ifstream = dynamic_cast<std::ifstream*>(&input);
      assert(ifstream != NULL);

      Imf::StdIFStream iifstream(*ifstream, "unknown file name");

      Imf::RgbaInputFile file(iifstream);
			Imath::Box2i dw = file.dataWindow();

			W = dw.max.x - dw.min.x + 1;
			H = dw.max.y - dw.min.y + 1;

			Imf::Array2D<Imf::Rgba> pixels;
			pixels.resizeErase(H, W);

			file.setFrameBuffer (&pixels[0][0] - dw.min.x - dw.min.y * W, 1, W);
			file.readPixels (dw.min.y, dw.max.y);

			pix = new FType[W*H*3] ;

			switch(nC)
			{
				case 3: for(int i=0;i<H;++i)
							  for(int j=0;j<W;++j)
							  {
								  pix[3*(j+i*W)+0] = pixels[H-i-1][j].r ;
								  pix[3*(j+i*W)+1] = pixels[H-i-1][j].g ;
								  pix[3*(j+i*W)+2] = pixels[H-i-1][j].b ;
							  }
						  break ;

				case 1: for(int i=0;i<H;++i)
							  for(int j=0;j<W;++j)
								  pix[j+i*W] = 0.3*pixels[H-i-1][j].r + 0.59*pixels[H-i-1][j].g + 0.11*pixels[H-i-1][j].b ;
						  break ;
				default:
						  throw std::runtime_error("Unexpected case in EXR_IO.") ;
			}
			return true ;
		}

		static bool SaveEXR(const char *filename,int W,int H, const FType *pix,int nC=3)
		{
			Imf::Array2D<Imf::Rgba> pixels(H,W);

			/* Convert separated channel representation to per pixel representation */

			switch(nC)
			{
				case 3:
					for (int row=0;row<H;row++) 
						for(int i=0;i<W;i++) 
						{
							Imf::Rgba &p = pixels[H-row-1][i];

							p.r = pix[3*(i+row*W)+0] ;
							p.g = pix[3*(i+row*W)+1] ;
							p.b = pix[3*(i+row*W)+2] ;
							p.a = 1.0 ;
						}
					break ;

				case 1:
					for (int row=0;row<H;row++) 
						for(int i=0;i<W;i++) 
						{
							Imf::Rgba &p = pixels[H-row-1][i];

							p.r = pix[i+row*W] ;
							p.g = pix[i+row*W] ;
							p.b = pix[i+row*W] ;
							p.a = FType(1.0) ;
						}
					break ;
				default:
					throw std::runtime_error("Unexpected case in EXR_IO.") ;
			}

			Imf::RgbaOutputFile file(filename, W, H, Imf::WRITE_RGBA);
			file.setFrameBuffer(&pixels[0][0], 1, W);
			file.writePixels(H);

			return true ;
		}
};

typedef t_EXR_IO<float> EXR_IO ;


#endif
