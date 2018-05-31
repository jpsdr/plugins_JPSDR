#include "resample.h"
#include "AutoYUY2.h"
#include "nnedi3.h"
#include "aWarpSharp.h"

ThreadPoolInterface *poolInterface;

bool aWarpSharp_Enable_SSE2,aWarpSharp_Enable_SSE41,aWarpSharp_Enable_AVX;

const AVS_Linkage *AVS_linkage = nullptr;


#define PLUGINS_JPSDR_VERSION "Plugins JPSDR 2.2.0"

/*
  threshold : int, default value : 4
     Value for threshold between odd/even lines to consider pixel interlaced.
  mode : int, default value -1.
     Conversion mode YV12 to YUY2.
	 -1 : Automatic detection.
	 0 : Progressive.
	 1 : Interlaced.
	 2 : Test mode (put white/colored pixel on part detected interlaced).
  output : int, default value 1
     0 : Output YUY2
	 1 : Output YV16
  threads : int, default value 0
     0 : Automatic.
	 Other value specify the number of threads, max is MAX_MT_THREADS
*/
AVSValue __cdecl Create_AutoYUY2(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("AutoYUY2: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if (!vi.IsYV12())
		env->ThrowError("AutoYUY2: Input format must be YV12 or I420");
	if ((vi.width & 1)!=0)
		env->ThrowError("AutoYUY2: Input width must be a multiple of 2.");
	if ((vi.height & 3)!=0)
		env->ThrowError("AutoYUY2: Input height must be a multiple of 4.");
	if (vi.height < 8)
		env->ThrowError("AutoYUY2: Input height must be at least 8.");

	const int thrs=args[1].AsInt(4);
	const int mode=args[2].AsInt(-1);
	const int output=args[3].AsInt(1);
	const int threads=args[4].AsInt(0);
	const bool LogicalCores=args[5].AsBool(true);
	const bool MaxPhysCores=args[6].AsBool(true);
	const bool SetAffinity=args[7].AsBool(false);
	const bool sleep = args[8].AsBool(false);
	int prefetch=args[9].AsInt(0);

	if ((mode<-1) || (mode>2))
		env->ThrowError("AutoYUY2: [mode] must be -1 (Automatic), 0 (Progessive) , 1 (Interlaced) or 2 (Test).");
	if ((output<0) || (output>1))
		env->ThrowError("AutoYUY2: [output] must be 0 (YUY2) or 1 (YV16)");
	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("AutoYUY2: [threads] must be between 0 and %ld.",MAX_MT_THREADS);

	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL)) env->ThrowError("AutoYUY2: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	uint8_t threads_number=1;

	if (threads!=1)
	{
		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("AutoYUY2: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("AutoYUY2: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("AutoYUY2: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("AutoYUY2: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("AutoYUY2: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new AutoYUY2(args[0].AsClip(), thrs, mode, output, threads_number, sleep, env);
}


AVSValue __cdecl Create_nnedi3(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip())
		env->ThrowError("nnedi3: arg 0 must be a clip!");
	VideoInfo vi = args[0].AsClip()->GetVideoInfo();
	
  const bool avsp=env->FunctionExists("ConvertBits");
  const bool RGB32=vi.IsRGB32();
  const bool RGB48=vi.IsRGB48();
  const bool RGB64=vi.IsRGB64();

	if (avsp)
	{
		if (!vi.IsPlanar() && !vi.IsYUY2() && !vi.IsRGB24() && !RGB32 && !RGB48 && !RGB64)
			env->ThrowError("nnedi3: only planar, YUY2, RGB64, RGB48, RGB32 and RGB24 input are supported!");				
	}
	else
	{
		if (!vi.IsPlanar() && !vi.IsYUY2() && !vi.IsRGB24())
			env->ThrowError("nnedi3: only planar, YUY2 and RGB24 input are supported!");				
	}

	const int threads=args[11].AsInt(0);
	const bool LogicalCores=args[14].AsBool(true);
	const bool MaxPhysCores=args[15].AsBool(true);
	const bool SetAffinity=args[16].AsBool(false);
	const bool sleep = args[18].AsBool(false);
	int prefetch = args[19].AsInt(0);

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("nnedi3: [threads] must be between 0 and %ld.",MAX_MT_THREADS);

	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL)) env->ThrowError("nnedi3: [prefetch] must be between 0 and %d.", MAX_THREAD_POOL);
			
	const bool dh = args[2].AsBool(false);
	if (((vi.height&1)!=0) && !dh)
		env->ThrowError("nnedi3: height must be mod 2 when dh=false (%d)!", vi.height);

	uint8_t threads_number=1;

	if (threads!=1)
	{
		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("nnedi3: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("nnedi3: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("nnedi3: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("nnedi3: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("nnedi3: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	if (!vi.IsY8())
	{
		if (RGB32 || RGB48 || RGB64)
		{
			AVSValue v=args[0].AsClip();
			
			if (RGB32 || RGB64) v=env->Invoke("ConvertToPlanarRGBA",v).AsClip();
			else v=env->Invoke("ConvertToPlanarRGB",v).AsClip();
			v= new nnedi3(v.AsClip(),args[1].AsInt(-1),args[2].AsBool(false),
				args[3].AsBool(true),args[4].AsBool(true),args[5].AsBool(true),args[17].AsBool(true),
				args[6].AsInt(6),args[7].AsInt(1),args[8].AsInt(1),args[9].AsInt(0),args[10].AsInt(2),
				threads_number,args[12].AsInt(0),args[13].AsInt(15),args[18].AsBool(false),args[20].AsInt(1),avsp,env);
			if (RGB32) return env->Invoke("ConvertToRGB32",v).AsClip();
			else
			{
				if (RGB48) return env->Invoke("ConvertToRGB48",v).AsClip();
				else return env->Invoke("ConvertToRGB64",v).AsClip();
			}
		}
		else return new nnedi3(args[0].AsClip(),args[1].AsInt(-1),args[2].AsBool(false),
				args[3].AsBool(true),args[4].AsBool(true),args[5].AsBool(true),args[17].AsBool(true),
				args[6].AsInt(6),args[7].AsInt(1),args[8].AsInt(1),args[9].AsInt(0),args[10].AsInt(2),
				threads_number,args[12].AsInt(0),args[13].AsInt(15),args[18].AsBool(false),args[20].AsInt(1),avsp,env);			
	}
	else
		return new nnedi3(args[0].AsClip(),args[1].AsInt(-1),args[2].AsBool(false),
				args[3].AsBool(true),false,false,false,args[6].AsInt(6),args[7].AsInt(1),args[8].AsInt(1),
				args[9].AsInt(0),args[10].AsInt(2),threads_number,args[12].AsInt(0),args[13].AsInt(15),args[18].AsBool(false),
				args[20].AsInt(1),avsp,env);
}


AVSValue __cdecl Create_nnedi3_rpow2(AVSValue args, void* user_data, IScriptEnvironment *env)
{
	if (!args[0].IsClip()) env->ThrowError("nnedi3_rpow2: arg 0 must be a clip!");
	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

  const bool avsp=env->FunctionExists("ConvertBits");  
  const bool RGB32=vi.IsRGB32();
  const bool RGB48=vi.IsRGB48();
  const bool RGB64=vi.IsRGB64();
  const bool RGB24=vi.IsRGB24();  
  const bool isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
  const bool grey = vi.IsY();  
  const bool isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();
	
	if (avsp)
	{
		if (!vi.IsPlanar() && !vi.IsYUY2() && !RGB24 && !RGB32 && !RGB48 && !RGB64)
			env->ThrowError("nnedi3: only planar, YUY2, RGB32, RGB48, RGB64 and RGB24 input are supported!");				
	}
	else
	{
		if (!vi.IsPlanar() && !vi.IsYUY2() && !RGB24)
			env->ThrowError("nnedi3: only planar, YUY2 and RGB24 input are supported!");				
	}
  
	if ((vi.IsYUY2() || vi.IsYV16()|| vi.IsYV12() || vi.IsYV411()) && (vi.width&3))
		env->ThrowError("nnedi3_rpow2: for YV12, YV411, YUY2 and YV16 input width must be mod 4 (%d)!", vi.width);
	const int rfactor = args[1].AsInt(-1);
	const int nsize = args[2].AsInt(0);
	const int nns = args[3].AsInt(3);
	const int qual = args[4].AsInt(1);
	const int etype = args[5].AsInt(0);
	const int pscrn = args[6].AsInt(2);
	const char *cshift = args[7].AsString("");
	const int fwidth = args[8].IsInt() ? args[8].AsInt() : rfactor*vi.width;
	const int fheight = args[9].IsInt() ? args[9].AsInt() : rfactor*vi.height;
	const float ep0 = (float)(args[10].IsFloat() ? args[10].AsFloat() : -FLT_MAX);
	const float ep1 = (float)(args[11].IsFloat() ? args[11].AsFloat() : -FLT_MAX);
	const int threads = args[12].AsInt(0);
	const int opt = args[13].AsInt(0);
	const int fapprox = args[14].AsInt(15);
	const bool chroma_shift_resize = args[15].AsBool(true);
	const bool mpeg2_chroma = args[16].AsBool(true);
	const bool LogicalCores = args[17].AsBool(true);
	const bool MaxPhysCores = args[18].AsBool(true);
	const bool SetAffinity = args[19].AsBool(false);
	const int threads_rs = args[20].AsInt(0);
	const bool LogicalCores_rs = args[21].AsBool(true);
	const bool MaxPhysCores_rs = args[22].AsBool(true);
	const bool SetAffinity_rs = args[23].AsBool(false);
	const bool sleep = args[24].AsBool(false);
	int prefetch = args[25].AsInt(0);
	int range_mode = args[26].AsInt(1);

	if ((rfactor<2) || (rfactor>1024)) env->ThrowError("nnedi3_rpow2: 2 <= rfactor <= 1024, and rfactor be a power of 2!\n");
	int rf=1,ct=0;

	while (rf<rfactor)
	{
		rf*=2;
		ct++;
	}

	if (rf!=rfactor)
		env->ThrowError("nnedi3_rpow2: 2 <= rfactor <= 1024, and rfactor be a power of 2!\n");
	if (nsize < 0 || nsize >= NUM_NSIZE)
		env->ThrowError("nnedi3_rpow2: nsize must be in [0,%d]!\n", NUM_NSIZE-1);
	if (nns < 0 || nns >= NUM_NNS)
		env->ThrowError("nnedi3_rpow2: nns must be in [0,%d]!\n", NUM_NNS-1);
	if (qual < 1 || qual > 2)
		env->ThrowError("nnedi3_rpow2: qual must be set to 1 or 2!\n");
	if (threads < 0 || threads > MAX_MT_THREADS)
		env->ThrowError("nnedi3_rpow2: 0 <= threads <= %d!\n",MAX_MT_THREADS);
	if (threads_rs < 0 || threads_rs > MAX_MT_THREADS)
		env->ThrowError("nnedi3_rpow2: 0 <= threads_rs <= %d!\n",MAX_MT_THREADS);
	if (opt < 0 || opt > 7)
		env->ThrowError("nnedi3_rpow2: opt must be in [0,7]!\n");
	if (fapprox < 0 || fapprox > 15)
		env->ThrowError("nnedi3_rpow2: fapprox must be [0,15]!\n");

	if ((range_mode<0) || (range_mode>4)) env->ThrowError("nnedi3_rpow2: [range] must be between 0 and 4!");

	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL)) env->ThrowError("nnedi3_rpow2: [prefetch] must be between 0 and %d.", MAX_THREAD_POOL);

	uint8_t threads_number=1;

	if (threads!=1)
	{
		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("nnedi3_rpow2: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("nnedi3_rpow2: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("nnedi3_rpow2: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("nnedi3_rpow2: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("nnedi3_rpow2: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	AVSValue v = args[0].AsClip();

	const bool FTurnL=(!avsp) && (env->FunctionExists("FTurnLeft") && ((env->GetCPUFlags() & CPUF_SSE2)!=0));
	const bool FTurnR=(!avsp) && (env->FunctionExists("FTurnRight") && ((env->GetCPUFlags() & CPUF_SSE2)!=0));
	const bool SplineMT=env->FunctionExists("Spline36ResizeMT");

	auto turnRightFunction = (FTurnR) ? "FTurnRight" : "TurnRight";
	auto turnLeftFunction =  (FTurnL) ? "FTurnLeft" : "TurnLeft";
	auto Spline36 = (SplineMT) ? "Spline36ResizeMT" : "Spline36Resize";

	uint8_t plane_range[PLANE_MAX];

	if ((range_mode!=1) && (range_mode!=4))
	{
		if (vi.IsYUV())
		{
			plane_range[0]=2;
			plane_range[1]=3;
			plane_range[2]=3;
		}
		else
		{
			if (grey)
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=(range_mode==0) ? 2 : range_mode;
			}
			else
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=1;
			}
		}
	}
	else
	{
		if (vi.IsRGB()) range_mode=1;

		for (uint8_t i=0; i<3; i++)
			plane_range[i]=range_mode;
	}
	plane_range[3]=1;
	range_mode=plane_range[0];

	try 
	{
		double Y_hshift=0.0,Y_vshift=0.0,C_hshift=0.0,C_vshift=0.0;

		const bool do_resize=(cshift[0]!=0) || vi.Is420();

		AVSValue vv,vu,va;

		if (RGB24 || vi.Is444() || vi.IsY() || RGB32 || RGB48 || RGB64 || isRGBPfamily)
		{
			Y_hshift = Y_vshift = -0.5;
		}
		else
		{
			Y_hshift = -0.5*(rf-1);
			Y_vshift = -0.5;

			C_hshift=Y_hshift;
			C_vshift=Y_vshift;

			if (vi.Is420())
			{
				// Correct chroma shift (it's always 1/2 pixel upwards).
				C_vshift-=0.5;

				C_vshift/=2.0;
				C_hshift/=2.0;

				C_hshift-=0.25*(rf-1);

				// Correct resize chroma position if YV12 has MPEG2 chroma subsampling
				if (chroma_shift_resize && mpeg2_chroma && (fwidth!=vi.width))
					C_hshift+=0.25*rf*(1.0-(double)vi.width/(double)fwidth);
			}
			else
			{
				if (vi.IsYV411())
				{
					C_hshift/=4.0;
					C_hshift-=0.375*(rf-1);

				// Correct resize chroma position
				if (chroma_shift_resize && (fwidth!=vi.width))
					C_hshift+=0.375*rf*(1.0-(double)vi.width/(double)fwidth);
				}
				else
				{
					C_hshift/=2.0;
					C_hshift-=0.25*(rf-1);

					//YV16 always has MPEG2 chroma subsampling
					if (chroma_shift_resize && (fwidth!=vi.width))
						C_hshift+=0.25*rf*(1.0-(double)vi.width/(double)fwidth);
				}
			}
		}

		if (RGB24 || vi.Is444() || vi.IsY() || RGB32 || RGB48 || RGB64 || isRGBPfamily)
		{
			if (RGB24 || RGB48)
			{
				if (avsp) v=env->Invoke("ConvertToPlanarRGB",v).AsClip();
				else
				{
					AVSValue sargs[3] = {v,"Y8",0};
					
					vu=env->Invoke("ShowRed",AVSValue(sargs,2)).AsClip();
					vv=env->Invoke("ShowGreen",AVSValue(sargs,2)).AsClip();
					v=env->Invoke("ShowBlue",AVSValue(sargs,2)).AsClip();
					sargs[0]=vu; sargs[1]=vv; sargs[2]=v;
					v=env->Invoke("Interleave",AVSValue(sargs,3)).AsClip();
				}					
			}
			
			if (RGB32 || RGB64) v=env->Invoke("ConvertToPlanarRGBA",v).AsClip();

			const bool UV_process=!(vi.IsY() || (RGB24 && !avsp));

			const int range_=(do_resize) ? 1 : range_mode;

			for (int i=0; i<ct; i++)
			{
				v = env->Invoke(turnRightFunction,v).AsClip();
				v = new nnedi3(v.AsClip(),i==0?1:0,true,true,UV_process,UV_process,isAlphaChannel || RGB32,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,1,avsp,env);
				v = env->Invoke(turnLeftFunction,v).AsClip();
				v = new nnedi3(v.AsClip(),i==0?1:0,true,true,UV_process,UV_process,isAlphaChannel || RGB32,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,(i==(ct-1))?range_:1,avsp,env);
			}
		}
		else
		{
			if (avsp && !vi.IsYUY2())
			{
				AVSValue sargs[2] = {v,"U"};
				
				vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
				sargs[1]="V";				
				vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
				if (isAlphaChannel)
				{
					sargs[1]="A";
					va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();					
				}
				sargs[1]="Y";
				v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();				
			}
			else
			{
				vu = env->Invoke("UtoY8",v).AsClip();
				vv = env->Invoke("VtoY8",v).AsClip();
				v = env->Invoke("ConvertToY8",v).AsClip();				
			}

			int range_=(do_resize) ? 1 : plane_range[0];

			for (int i=0; i<ct; i++)
			{
				v = env->Invoke(turnRightFunction,v).AsClip();
				// always use field=1 to keep chroma/luma horizontal alignment
				v = new nnedi3(v.AsClip(),1,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,1,avsp,env);
				v = env->Invoke(turnLeftFunction,v).AsClip();
				v = new nnedi3(v.AsClip(),i==0?1:0,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,(i==(ct-1))?range_:1,avsp,env);
			}

			range_=(do_resize) ? 1 : plane_range[1];
			for (int i=0; i<ct; i++)
			{
				vu = env->Invoke(turnRightFunction,vu).AsClip();
				// always use field=1 to keep chroma/luma horizontal alignment
				vu = new nnedi3(vu.AsClip(),1,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,1,avsp,env);
				vu = env->Invoke(turnLeftFunction,vu).AsClip();
				vu = new nnedi3(vu.AsClip(),i==0?1:0,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,(i==(ct-1))?range_:1,avsp,env);
			}

			range_=(do_resize) ? 1 : plane_range[2];
			for (int i=0; i<ct; i++)
			{
				vv = env->Invoke(turnRightFunction,vv).AsClip();
				// always use field=1 to keep chroma/luma horizontal alignment
				vv = new nnedi3(vv.AsClip(),1,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,1,avsp,env);
				vv = env->Invoke(turnLeftFunction,vv).AsClip();
				vv = new nnedi3(vv.AsClip(),i==0?1:0,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,(i==(ct-1))?range_:1,avsp,env);
			}

			range_=(do_resize) ? 1 : plane_range[3];
			if (isAlphaChannel)
			{
				for (int i=0; i<ct; i++)
				{
					va = env->Invoke(turnRightFunction,va).AsClip();
					// always use field=1 to keep chroma/luma horizontal alignment
					va = new nnedi3(va.AsClip(),1,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,1,avsp,env);
					va = env->Invoke(turnLeftFunction,va).AsClip();
					va = new nnedi3(va.AsClip(),i==0?1:0,true,true,false,false,false,nsize,nns,qual,etype,pscrn,threads_number,opt,fapprox,sleep,(i==(ct-1))?range_:1,avsp,env);
				}				
			}
		}

		if (cshift[0]!=0)
		{
			const bool use_rs_mt=((_strnicmp(cshift,"pointresizemt",13)==0) || (_strnicmp(cshift,"bilinearresizemt",16)==0)
				|| (_strnicmp(cshift,"bicubicresizemt",15)==0) || (_strnicmp(cshift,"lanczosresizemt",15)==0)
				|| (_strnicmp(cshift,"lanczos4resizemt",16)==0) || (_strnicmp(cshift,"blackmanresizemt",16)==0)
				|| (_strnicmp(cshift,"spline16resizemt",16)==0) || (_strnicmp(cshift,"spline36resizemt",16)==0)
				|| (_strnicmp(cshift,"spline64resizemt",16)==0) || (_strnicmp(cshift,"gaussresizemt",13)==0)
				|| (_strnicmp(cshift,"sincresizemt",12)==0));

			int type = 0;
			
			if ((_strnicmp(cshift,"blackmanresize",14)==0) || (_strnicmp(cshift,"lanczosresize",13)==0)
				|| (_strnicmp(cshift,"sincresize",10)==0)) type=1;
			else
			{
				if (_strnicmp(cshift,"gaussresize",11)==0) type=2;
				else
				{
					if (_strnicmp(cshift,"bicubicresize",13)==0) type=3;
				}
			}
			
			if ((type==0) || ((type!=3) && (ep0==-FLT_MAX)) ||
				((type==3) && (ep0==-FLT_MAX) && (ep1==-FLT_MAX)))
			{
				AVSValue sargs[14] = { v, fwidth, fheight, Y_hshift, Y_vshift, 
					vi.width*rfactor, vi.height*rfactor,threads_rs,LogicalCores_rs,MaxPhysCores_rs,SetAffinity_rs,sleep,
					prefetch,range_mode };
				const char *nargs[14] = { 0, 0, 0, "src_left", "src_top", 
					"src_width", "src_height","threads","logicalCores","MaxPhysCore","SetAffinity","sleep","prefetch","range" };
				const uint8_t nbarg=(use_rs_mt) ? 14:7;

				v=env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();

				if (!(RGB24 || vi.Is444() || vi.IsY() || RGB32 || RGB48 || RGB64 || isRGBPfamily))
				{
					if (isAlphaChannel)
					{
						sargs[0]=va;
						sargs[13]=plane_range[3];
						va=env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();
					}
					
					sargs[3]=C_hshift;
					sargs[4]=C_vshift;

					if (vi.Is420())
					{
						sargs[1]=fwidth >> 1;
						sargs[2]=fheight >> 1;
						sargs[5]=(vi.width*rfactor) >> 1;
						sargs[6]=(vi.height*rfactor) >> 1;
					}
					else
					{
						if (vi.IsYV411())
						{
							sargs[1]=fwidth >> 2;
							sargs[5]=(vi.width*rfactor) >> 2;
						}
						else
						{
							sargs[1]=fwidth >> 1;
							sargs[5]=(vi.width*rfactor) >> 1;
						}
					}

					sargs[0]=vu;
					sargs[13]=plane_range[1];
					vu = env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();
					sargs[0]=vv;
					sargs[13]=plane_range[2];
					vv = env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();

					AVSValue ytouvargs[4] = {vu,vv,v,va};
					if (isAlphaChannel) v=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
					else v=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();

					if (vi.IsYUY2()) v=env->Invoke("ConvertToYUY2",v).AsClip();
				}
				else
				{
					if (RGB24)
					{
						if (avsp) v=env->Invoke("ConvertToRGB24",v).AsClip();
						else
						{
							sargs[0]=v; sargs[1]=3;
							sargs[2]=0;
							vu=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[2]=1;
							vv=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[2]=2;
							v=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[0]=vu; sargs[1]=vv; sargs[2]=v; sargs[3]="RGB24";
							v=env->Invoke("MergeRGB",AVSValue(sargs,4)).AsClip();							
						}
					}
					if (RGB32) v=env->Invoke("ConvertToRGB32",v).AsClip();
					if (RGB48) v=env->Invoke("ConvertToRGB48",v).AsClip();
					if (RGB64) v=env->Invoke("ConvertToRGB64",v).AsClip();
				}
			}
			else if ((type!=3) || (min(ep0,ep1)==-FLT_MAX))
			{
				AVSValue sargs[15] = { v, fwidth, fheight, Y_hshift, Y_vshift, 
					vi.width*rfactor, vi.height*rfactor, type==1?AVSValue((int)(ep0+0.5f)):
					(type==2?ep0:max(ep0,ep1)),threads_rs,LogicalCores_rs,MaxPhysCores_rs,SetAffinity_rs,sleep,prefetch,range_mode };
				const char *nargs[15] = { 0, 0, 0, "src_left", "src_top", 
					"src_width", "src_height", type==1?"taps":(type==2?"p":(max(ep0,ep1)==ep0?"b":"c")),
					"threads","logicalCores","MaxPhysCore","SetAffinity","sleep","prefetch","range" };
				const uint8_t nbarg=(use_rs_mt) ? 15:8;

				v=env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();

				if (!(RGB24 || vi.Is444() || vi.IsY() || RGB32 || RGB48 || RGB64 || isRGBPfamily))
				{
					if (isAlphaChannel)
					{
						sargs[0]=va;
						sargs[14]=plane_range[3];
						va=env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();
					}
					
					sargs[3]=C_hshift;
					sargs[4]=C_vshift;

					if (vi.Is420())
					{
						sargs[1]=fwidth >> 1;
						sargs[2]=fheight >> 1;
						sargs[5]=(vi.width*rfactor) >> 1;
						sargs[6]=(vi.height*rfactor) >> 1;
					}
					else
					{
						if (vi.IsYV411())
						{
							sargs[1]=fwidth >> 2;
							sargs[5]=(vi.width*rfactor) >> 2;
						}
						else
						{
							sargs[1]=fwidth >> 1;
							sargs[5]=(vi.width*rfactor) >> 1;
						}
					}

					sargs[0]=vu;
					sargs[14]=plane_range[1];
					vu = env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();
					sargs[0]=vv;
					sargs[14]=plane_range[2];
					vv = env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();

					AVSValue ytouvargs[4] = {vu,vv,v,va};
					if (isAlphaChannel) v=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
					else v=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();

					if (vi.IsYUY2()) v=env->Invoke("ConvertToYUY2",v).AsClip();
				}
				else
				{
					if (RGB24)
					{
						if (avsp) v=env->Invoke("ConvertToRGB24",v).AsClip();
						else
						{
							sargs[0]=v; sargs[1]=3;
							sargs[2]=0;
							vu=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[2]=1;
							vv=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[2]=2;
							v=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[0]=vu; sargs[1]=vv; sargs[2]=v; sargs[3]="RGB24";
							v=env->Invoke("MergeRGB",AVSValue(sargs,4)).AsClip();							
						}
					}
					if (RGB32) v=env->Invoke("ConvertToRGB32",v).AsClip();
					if (RGB48) v=env->Invoke("ConvertToRGB48",v).AsClip();
					if (RGB64) v=env->Invoke("ConvertToRGB64",v).AsClip();
				}
			}
			else
			{
				AVSValue sargs[16] = { v, fwidth, fheight, Y_hshift, Y_vshift, 
					vi.width*rfactor, vi.height*rfactor, ep0, ep1,threads_rs,LogicalCores_rs,MaxPhysCores_rs,
					SetAffinity_rs,sleep,prefetch,range_mode };
				const char *nargs[16] = { 0, 0, 0, "src_left", "src_top", 
					"src_width", "src_height", "b", "c", "threads","logicalCores","MaxPhysCore","SetAffinity",
					"sleep","prefetch","range" };
				const uint8_t nbarg=(use_rs_mt) ? 16:9;

				v = env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();

				if (!(RGB24 || vi.Is444() || vi.IsY() || RGB32 || RGB48 || RGB64 || isRGBPfamily))
				{
					if (isAlphaChannel)
					{
						sargs[0]=va;
						sargs[15]=plane_range[3];
						va=env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();
					}
					
					sargs[3]=C_hshift;
					sargs[4]=C_vshift;

					if (vi.Is420())
					{
						sargs[1]=fwidth >> 1;
						sargs[2]=fheight >> 1;
						sargs[5]=(vi.width*rfactor) >> 1;
						sargs[6]=(vi.height*rfactor) >> 1;
					}
					else
					{
						if (vi.IsYV411())
						{
							sargs[1]=fwidth >> 2;
							sargs[5]=(vi.width*rfactor) >> 2;
						}
						else
						{
							sargs[1]=fwidth >> 1;
							sargs[5]=(vi.width*rfactor) >> 1;
						}
					}

					sargs[0]=vu;
					sargs[15]=plane_range[1];
					vu = env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();
					sargs[0]=vv;
					sargs[15]=plane_range[2];
					vv = env->Invoke(cshift,AVSValue(sargs,nbarg),nargs).AsClip();

					AVSValue ytouvargs[4] = {vu,vv,v,va};
					if (isAlphaChannel) v=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
					else v=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();

					if (vi.IsYUY2()) v=env->Invoke("ConvertToYUY2",v).AsClip();
				}
				else
				{
					if (RGB24)
					{
						if (avsp) v=env->Invoke("ConvertToRGB24",v).AsClip();
						else
						{
							sargs[0]=v; sargs[1]=3;
							sargs[2]=0;
							vu=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[2]=1;
							vv=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[2]=2;
							v=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
							sargs[0]=vu; sargs[1]=vv; sargs[2]=v; sargs[3]="RGB24";
							v=env->Invoke("MergeRGB",AVSValue(sargs,4)).AsClip();							
						}
					}
					if (RGB32) v=env->Invoke("ConvertToRGB32",v).AsClip();
					if (RGB48) v=env->Invoke("ConvertToRGB48",v).AsClip();
					if (RGB64) v=env->Invoke("ConvertToRGB64",v).AsClip();
				}
			}
		}
		else
		{
			if (!(RGB24 || vi.Is444() || vi.IsY() || RGB32 || RGB48 || RGB64 || isRGBPfamily))
			{
				if (vi.Is420())
				{
					AVSValue sargs[14]={vu,(vi.width*rfactor)>>1,(vi.height*rfactor)>>1,0.0,-0.25,
						(vi.width*rfactor)>>1,(vi.height*rfactor)>>1,threads_rs,LogicalCores_rs,MaxPhysCores_rs,
						SetAffinity_rs,sleep,prefetch,plane_range[1]};
					const char *nargs[14]={0,0,0,"src_left","src_top","src_width","src_height","threads",
					"logicalCores","MaxPhysCore","SetAffinity","sleep","prefetch","range" };
					const uint8_t nbarg=(SplineMT) ? 14:7;

					vu = env->Invoke(Spline36,AVSValue(sargs,nbarg),nargs).AsClip();
					sargs[0]=vv;
					sargs[13]=plane_range[2];
					vv = env->Invoke(Spline36,AVSValue(sargs,nbarg),nargs).AsClip();
				}

				AVSValue ytouvargs[4] = {vu,vv,v,va};
				if (isAlphaChannel) v=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
				else v=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();

				if (vi.IsYUY2()) v=env->Invoke("ConvertToYUY2",v).AsClip();
			}
			else
			{
				if (RGB24)
				{
					if (avsp) v=env->Invoke("ConvertToRGB24",v).AsClip();
					else
					{
						AVSValue sargs[4];

						sargs[0]=v; sargs[1]=3;
						sargs[2]=0;
						vu=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
						sargs[2]=1;
						vv=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
						sargs[2]=2;
						v=env->Invoke("SelectEvery",AVSValue(sargs,3)).AsClip();
						sargs[0]=vu; sargs[1]=vv; sargs[2]=v; sargs[3]="RGB24";
						v=env->Invoke("MergeRGB",AVSValue(sargs,4)).AsClip();						
					}
				}				
				if (RGB32) v=env->Invoke("ConvertToRGB32",v).AsClip();
				if (RGB48) v=env->Invoke("ConvertToRGB48",v).AsClip();
				if (RGB64) v=env->Invoke("ConvertToRGB64",v).AsClip();
			}
		}
	}
	catch (IScriptEnvironment::NotFound)
	{
		env->ThrowError("nnedi3_rpow2: error using env->invoke (function not found)!\n");
	}
	return v;
}


static bool is_cplace_mpeg2(const AVSValue &args, int pos)
{
  const char *cplace_0=args[pos].AsString("");
  const bool cplace_mpeg2_flag=(_stricmp(cplace_0, "MPEG2")==0);
  return (cplace_mpeg2_flag);
}


// thresh: 0..255
// blur:   0..?
// type:   0..1
// depth:  -128..127
// chroma modes:
// 0 - zero
// 1 - don't care
// 2 - copy
// 3 - process
// 4 - guide by luma - warp only
// remap from MarcFD's aWarpSharp: thresh=_thresh*256, blur=_blurlevel, type= (bm=0)->1, (bm=2)->0, depth=_depth*_blurlevel/2, chroma= 0->2, 1->4, 2->3
AVSValue __cdecl Create_aWarpSharp(AVSValue args, void *user_data, IScriptEnvironment *env)
{
	int threads,prefetch;
	bool LogicalCores,MaxPhysCores,SetAffinity,sleep;

	uint8_t threads_number=1;

	if (!args[0].IsClip()) env->ThrowError("aWarpSharpMT: arg 0 must be a clip!");
	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	const bool avsp=env->FunctionExists("ConvertBits");

	int thresh,blur,blurt,depth,depthC,blurV,depthV,depthVC,blurC,blurVC,threshC;

  switch ((int)(size_t)user_data)
  {
  case 0 :
	  {
	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aWarpSharp2: SSE2 capable CPU is required");

	  threads=args[14].AsInt(0);
	  LogicalCores=args[15].AsBool(true);
	  MaxPhysCores=args[16].AsBool(true);
	  SetAffinity=args[17].AsBool(false);
	  sleep = args[18].AsBool(false);
	  prefetch=args[19].AsInt(0);

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarpSharp2: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarpSharp2: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	  if (threads!=1)
	  {
		  if (!poolInterface->CreatePool(prefetch)) env->ThrowError("aWarpSharp2: Unable to create ThreadPool!");

		  threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		  if (threads_number==0) env->ThrowError("aWarpSharp2: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("aWarpSharp2: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("aWarpSharp2: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("aWarpSharp2: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	  }

	  thresh=args[1].AsInt(0x80);
	  blurt=args[3].AsInt(0);
	  args[2].Defined() ? blur=args[2].AsInt(-1) : blur=((blurt==0) ? 2:3);
	  depth=args[4].AsInt(16);
	  args[6].Defined() ? depthC=args[6].AsInt(128) : depthC=(vi.Is444() ? depth:(depth>>1));
	  args[8].Defined() ? blurV=args[8].AsInt(-1) : blurV=blur;
	  args[9].Defined() ? depthV=args[9].AsInt(128) : depthV=depth;
	  args[10].Defined() ? depthVC=args[10].AsInt(128) : depthVC=depthC;
	  args[11].Defined() ? blurC=args[11].AsInt(-1) : blurC=(blur+1)>>1;
	  args[12].Defined() ? blurVC=args[12].AsInt(-1) : blurVC=blurC;
	  args[13].Defined() ? threshC=args[13].AsInt(-1) : threshC=thresh;

    return new aWarpSharp(args[0].AsClip(),thresh,blur,blurt,depth,args[5].AsInt(4),depthC,is_cplace_mpeg2(args,7),
		blurV,depthV,depthVC,blurC,blurVC,threshC,threads_number,sleep,avsp,env);
	break;
	  }
  case 1 :
	  {
	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aWarpSharp: SSE2 capable CPU is required");

    blurt = (args[5].AsInt(2)!=2)?1:0;
    const int blurlevel = args[2].AsInt(2);
    const unsigned int cm = args[4].AsInt(1);
    static const char map[4] = {1,4,3,2};

	  threads=args[7].AsInt(0);
	  LogicalCores=args[8].AsBool(true);
	  MaxPhysCores=args[9].AsBool(true);
	  SetAffinity=args[10].AsBool(false);
	  sleep = args[11].AsBool(false);
	  prefetch=args[12].AsInt(0);

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarpSharp: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarpSharp: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	  if (threads!=1)
	  {
		  if (!poolInterface->CreatePool(prefetch)) env->ThrowError("aWarpSharp: Unable to create ThreadPool!");

		  threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		  if (threads_number==0) env->ThrowError("aWarpSharp: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("aWarpSharp: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("aWarpSharp: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("aWarpSharp: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	  }

	  thresh=int(args[3].AsFloat(0.5)*256.0);
	  blur=(blurt==1)?(blurlevel*3):blurlevel;
	  depth=int(args[1].AsFloat(16.0)*blurlevel*0.5);
	  depthC=vi.Is444() ? depth:(depth>>1);
	  blurV=blur;
	  depthV=depth;
	  depthVC=depthC;
	  blurC=(blur+1)>>1;
	  blurVC=blurC;
	  threshC=thresh;

    return new aWarpSharp(args[0].AsClip(),thresh,blur,blurt,depth,(cm<4)?map[cm]:-1,depthC,false,
		blurV,depthV,depthVC,blurC,blurVC,threshC,threads_number,sleep,avsp,env);
	break;
	  }
  case 2 :
	  {

	  threads=args[4].AsInt(0);
	  LogicalCores=args[5].AsBool(true);
	  MaxPhysCores=args[6].AsBool(true);
	  SetAffinity=args[7].AsBool(false);
	  sleep = args[8].AsBool(false);
	  prefetch=args[9].AsInt(0);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aSobel: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aSobel: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aSobel: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	  if (threads!=1)
	  {
		  if (!poolInterface->CreatePool(prefetch)) env->ThrowError("aSobel: Unable to create ThreadPool!");

		  threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		  if (threads_number==0) env->ThrowError("aSobel: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("aSobel: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("aSobel: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("aSobel: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	  }

	  thresh=args[1].AsInt(0x80);
	  args[3].Defined() ? threshC=args[3].AsInt(-1) : threshC=thresh;

	return new aSobel(args[0].AsClip(),thresh,args[2].AsInt(1),threshC,threads_number,sleep,avsp,env);
	break;
	  }
  case 3 :
	  {

	  threads=args[7].AsInt(0);
	  LogicalCores=args[8].AsBool(true);
	  MaxPhysCores=args[9].AsBool(true);
	  SetAffinity=args[10].AsBool(false);
	  sleep = args[11].AsBool(false);
	  prefetch=args[12].AsInt(0);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aBlur: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aBlur: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aBlur: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	  if (threads!=1)
	  {
		  if (!poolInterface->CreatePool(prefetch)) env->ThrowError("aBlur: Unable to create ThreadPool!");

		  threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		  if (threads_number==0) env->ThrowError("aBlur: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("aBlur: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("aBlur: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("aBlur: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	  }

	  blurt=args[2].AsInt(1);
	  args[1].Defined() ? blur=args[1].AsInt(-1) : blur=((blurt==0) ? 2:3);
	  args[4].Defined() ? blurV=args[4].AsInt(-1) : blurV=blur;
	  args[5].Defined() ? blurC=args[5].AsInt(-1) : blurC=(blur+1)>>1;
	  args[6].Defined() ? blurVC=args[6].AsInt(-1) : blurVC=blurC;

    return new aBlur(args[0].AsClip(),blur,blurt,args[3].AsInt(1),blurV,blurC,blurVC,threads_number,sleep,avsp,env);
	break;
	  }
  case 4 :
	  {

	  threads=args[8].AsInt(0);
	  LogicalCores=args[9].AsBool(true);
	  MaxPhysCores=args[10].AsBool(true);
	  SetAffinity=args[11].AsBool(false);
	  sleep = args[12].AsBool(false);
	  prefetch=args[13].AsInt(0);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aWarp: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarp: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarp: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	  if (threads!=1)
	  {
		  if (!poolInterface->CreatePool(prefetch)) env->ThrowError("aWarp: Unable to create ThreadPool!");

		  threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		  if (threads_number==0) env->ThrowError("aWarp: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("aWarp: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("aWarp: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("aWarp: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	  }

	  depth=args[2].AsInt(3);
	  args[4].Defined() ? depthC=args[4].AsInt(128) : depthC=(vi.Is444() ? depth:(depth>>1));
	  args[6].Defined() ? depthV=args[6].AsInt(128) : depthV=depth;
	  args[7].Defined() ? depthVC=args[7].AsInt(128) : depthVC=depthC;

    return new aWarp(args[0].AsClip(),args[1].AsClip(),depth,args[3].AsInt(4),depthC,is_cplace_mpeg2(args,5),
		depthV,depthVC,threads_number,sleep,avsp,env);
	break;
	  }
  case 5 :
	  {
	  threads=args[8].AsInt(0);
	  LogicalCores=args[9].AsBool(true);
	  MaxPhysCores=args[10].AsBool(true);
	  SetAffinity=args[11].AsBool(false);
	  sleep = args[12].AsBool(false);
	  prefetch=args[13].AsInt(0);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aWarp4: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarp4: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarp4: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	  if (threads!=1)
	  {
		  if (!poolInterface->CreatePool(prefetch)) env->ThrowError("aWarp4: Unable to create ThreadPool!");

		  threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		  if (threads_number==0) env->ThrowError("aWarp4: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,true,true,i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("aWarp4: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("aWarp4: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("aWarp4: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	  }

	  depth=args[2].AsInt(3);
	  args[4].Defined() ? depthC=args[4].AsInt(128) : depthC=(vi.Is444() ? depth:(depth>>1));
	  args[6].Defined() ? depthV=args[6].AsInt(128) : depthV=depth;
	  args[7].Defined() ? depthVC=args[7].AsInt(128) : depthVC=depthC;

    return new aWarp4(args[0].AsClip(),args[1].AsClip(),depth,args[3].AsInt(4),depthC,is_cplace_mpeg2(args,5),
		depthV,depthVC,threads_number,sleep,avsp,env);
	break;
	  }
  default : break;
  }
  return NULL;
}


extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
	AVS_linkage = vectors;

	poolInterface=ThreadPoolInterface::Init(0);

	if (!poolInterface->GetThreadPoolInterfaceStatus()) env->ThrowError("plugins_JPSDR: Error with the TheadPool status!");

	aWarpSharp_Enable_SSE2=(env->GetCPUFlags() & CPUF_SSE2)!=0;
	aWarpSharp_Enable_SSE41=(env->GetCPUFlags() & CPUF_SSE4_1)!=0;
	aWarpSharp_Enable_AVX=(env->GetCPUFlags() & CPUF_AVX)!=0;

	SetCPUMatrixClass((env->GetCPUFlags() & CPUF_SSE2)!=0,(env->GetCPUFlags() & CPUF_AVX)!=0,(env->GetCPUFlags() & CPUF_AVX2)!=0);

    env->AddFunction("AutoYUY2",
		"c[threshold]i[mode]i[output]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i",
		Create_AutoYUY2, 0);

	env->AddFunction("nnedi3", "c[field]i[dh]b[Y]b[U]b[V]b[nsize]i[nns]i[qual]i[etype]i[pscrn]i" \
		"[threads]i[opt]i[fapprox]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[A]b[sleep]b[prefetch]i[range]i", Create_nnedi3, 0);
	env->AddFunction("nnedi3_rpow2", "c[rfactor]i[nsize]i[nns]i[qual]i[etype]i[pscrn]i[cshift]s[fwidth]i" \
		"[fheight]i[ep0]f[ep1]f[threads]i[opt]i[fapprox]i[csresize]b[mpeg2]b[logicalCores]b[MaxPhysCore]b" \
		"[SetAffinity]b[threads_rs]i[logicalCores_rs]b[MaxPhysCore_rs]b[SetAffinity_rs]b[sleep]b[prefetch]i[range]i", Create_nnedi3_rpow2, 0);

	env->AddFunction("PointResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_PointResize, 0);
	env->AddFunction("BilinearResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_BilinearResize, 0);
	env->AddFunction("BicubicResizeMT", "c[target_width]i[target_height]i[b]f[c]f[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_BicubicResize, 0);
	env->AddFunction("LanczosResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_LanczosResize, 0);
	env->AddFunction("Lanczos4ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_Lanczos4Resize, 0);
	env->AddFunction("BlackmanResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_BlackmanResize, 0);
	env->AddFunction("Spline16ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_Spline16Resize, 0);
	env->AddFunction("Spline36ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_Spline36Resize, 0);
	env->AddFunction("Spline64ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_Spline64Resize, 0);
	env->AddFunction("GaussResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[p]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_GaussianResize, 0);
	env->AddFunction("SincResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i",
		FilteredResizeMT::Create_SincResize, 0);

	env->AddFunction("DeBilinearResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeBilinearResize, 0);
	env->AddFunction("DeBicubicResizeMT", "c[target_width]i[target_height]i[b]f[c]f[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeBicubicResize, 0);
	env->AddFunction("DeLanczosResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeLanczosResize, 0);
	env->AddFunction("DeLanczos4ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeLanczos4Resize, 0);
	env->AddFunction("DeBlackmanResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeBlackmanResize, 0);
	env->AddFunction("DeSpline16ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeSpline16Resize, 0);
	env->AddFunction("DeSpline36ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeSpline36Resize, 0);
	env->AddFunction("DeSpline64ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeSpline64Resize, 0);
	env->AddFunction("DeGaussResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[p]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeGaussianResize, 0);
	env->AddFunction("DeSincResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i",
		FilteredResizeMT::Create_DeSincResize, 0);

  env->AddFunction("aWarpSharp2", "c[thresh]i[blur]i[type]i[depth]i[chroma]i[depthC]i[cplace]s[blurV]i[depthV]i[depthVC]i" \
	  "[blurC]i[blurVC]i[threshC]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i", Create_aWarpSharp, (void*)0);
  env->AddFunction("aWarpSharp", "c[depth]f[blurlevel]i[thresh]f[cm]i[bm]i[show]b" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i", Create_aWarpSharp, (void*)1);
  env->AddFunction("aSobel", "c[thresh]i[chroma]i[threshC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i", Create_aWarpSharp, (void*)2);
  env->AddFunction("aBlur", "c[blur]i[type]i[chroma]i[blurV]i[blurC]i[blurVC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i", Create_aWarpSharp, (void*)3);
  env->AddFunction("aWarp", "cc[depth]i[chroma]i[depthC]i[cplace]s[depthV]i[depthVC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i", Create_aWarpSharp, (void*)4);
  env->AddFunction("aWarp4", "cc[depth]i[chroma]i[depthC]i[cplace]s[depthV]i[depthVC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i", Create_aWarpSharp, (void*)5);

	return "plugins JPSDR";	
}

