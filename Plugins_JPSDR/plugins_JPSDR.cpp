#include "./resample.h"
#include "./AutoYUY2.h"
#include "./nnedi3.h"
#include "./aWarpSharp.h"
#include "./HDRTools.h"

ThreadPoolInterface *poolInterface;

bool aWarpSharp_Enable_SSE2,aWarpSharp_Enable_SSE41,aWarpSharp_Enable_AVX;

const AVS_Linkage *AVS_linkage = nullptr;


#define PLUGINS_JPSDR_VERSION "Plugins JPSDR 3.2.3"

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

	const uint8_t pixelsize = (uint8_t)vi.ComponentSize();

	if ((pixelsize>16) || !vi.Is420())
		env->ThrowError("AutoYUY2: Input format must be YUV420, 8 to 16 bits");
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
	int thread_level=args[10].AsInt(6);

	if ((mode<-1) || (mode>2))
		env->ThrowError("AutoYUY2: [mode] must be -1 (Automatic), 0 (Progessive) , 1 (Interlaced) or 2 (Test).");
	if ((output<0) || (output>1))
		env->ThrowError("AutoYUY2: [output] must be 0 (YUY2) or 1 (YV16)");
	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("AutoYUY2: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("AutoYUY2: [ThreadLevel] must be between 1 and 7.");

	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL)) env->ThrowError("AutoYUY2: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

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
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("AutoYUY2: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("AutoYUY2: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("AutoYUY2: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new AutoYUY2(args[0].AsClip(), thrs, mode, output, threads_number,sleep, env);
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
	int prefetch = args[19].AsInt(0);
	int thread_level=args[21].AsInt(6);

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("nnedi3: [threads] must be between 0 and %ld.",MAX_MT_THREADS);

	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL)) env->ThrowError("nnedi3: [prefetch] must be between 0 and %d.", MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("nnedi3: [ThreadLevel] must be between 1 and 7.");
			
	const bool dh = args[2].AsBool(false);
	if (((vi.height&1)!=0) && !dh)
		env->ThrowError("nnedi3: height must be mod 2 when dh=false (%d)!", vi.height);

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

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
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("nnedi3: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("nnedi3: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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
	int thread_level=args[27].AsInt(6);
	int thread_level_rs=args[28].AsInt(6);

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
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("nnedi3_rpow2: [ThreadLevel] must be between 1 and 7.");
	if ((thread_level_rs<1) || (thread_level_rs>7))
		env->ThrowError("nnedi3_rpow2: [ThreadLevel_rs] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

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
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("nnedi3_rpow2: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("nnedi3_rpow2: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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
				AVSValue sargs[15] = { v, fwidth, fheight, Y_hshift, Y_vshift, 
					vi.width*rfactor, vi.height*rfactor,threads_rs,LogicalCores_rs,MaxPhysCores_rs,SetAffinity_rs,sleep,
					prefetch,range_mode,thread_level_rs };
				const char *nargs[15] = { 0, 0, 0, "src_left", "src_top", 
					"src_width", "src_height","threads","logicalCores","MaxPhysCore","SetAffinity","sleep",
					"prefetch","range","ThreadLevel" };
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
				AVSValue sargs[16] = { v, fwidth, fheight, Y_hshift, Y_vshift, 
					vi.width*rfactor, vi.height*rfactor, type==1?AVSValue((int)(ep0+0.5f)):
					(type==2?ep0:max(ep0,ep1)),threads_rs,LogicalCores_rs,MaxPhysCores_rs,SetAffinity_rs,
					sleep,prefetch,range_mode,thread_level_rs };
				const char *nargs[16] = { 0, 0, 0, "src_left", "src_top", 
					"src_width", "src_height", type==1?"taps":(type==2?"p":(max(ep0,ep1)==ep0?"b":"c")),
					"threads","logicalCores","MaxPhysCore","SetAffinity","sleep","prefetch","range","ThreadLevel" };
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
				AVSValue sargs[17] = { v, fwidth, fheight, Y_hshift, Y_vshift, 
					vi.width*rfactor, vi.height*rfactor, ep0, ep1,threads_rs,LogicalCores_rs,MaxPhysCores_rs,
					SetAffinity_rs,sleep,prefetch,range_mode,thread_level_rs };
				const char *nargs[17] = { 0, 0, 0, "src_left", "src_top", 
					"src_width", "src_height", "b", "c", "threads","logicalCores","MaxPhysCore","SetAffinity",
					"sleep","prefetch","range","ThreadLevel" };
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
					AVSValue sargs[15]={vu,(vi.width*rfactor)>>1,(vi.height*rfactor)>>1,0.0,-0.25,
						(vi.width*rfactor)>>1,(vi.height*rfactor)>>1,threads_rs,LogicalCores_rs,MaxPhysCores_rs,
						SetAffinity_rs,sleep,prefetch,plane_range[1],thread_level_rs};
					const char *nargs[15]={0,0,0,"src_left","src_top","src_width","src_height","threads",
					"logicalCores","MaxPhysCore","SetAffinity","sleep","prefetch","range","ThreadLevel" };
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
	int threads,prefetch,thread_level;
	bool LogicalCores,MaxPhysCores,SetAffinity,sleep;

	uint8_t threads_number=1;

	if (!args[0].IsClip()) env->ThrowError("aWarpSharpMT: arg 0 must be a clip!");
	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	const bool avsp=env->FunctionExists("ConvertBits");

	int thresh,blur,blurt,depth,depthC,blurV,depthV,depthVC,blurC,blurVC,threshC;

	const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
		BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

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
	  thread_level=args[20].AsInt(6);

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarpSharp2: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarpSharp2: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("aWarpSharp2: [ThreadLevel] must be between 1 and 7.");

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
						  if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,
							  MaxPhysCores,true,true,TabLevel[thread_level],i))
						  {
							  poolInterface->DeAllocateAllThreads(true);
							  env->ThrowError("aWarpSharp2: Error with the TheadPool while allocating threadpool!");
						  }
						  Offset+=delta;
					  }
				  }
				  else
				  {
					  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					  {
						  poolInterface->DeAllocateAllThreads(true);
						  env->ThrowError("aWarpSharp2: Error with the TheadPool while allocating threadpool!");
					  }
				  }
			  }
			  else
			  {
				  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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
	  thread_level=args[13].AsInt(6);

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarpSharp: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarpSharp: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("aWarpSharp: [ThreadLevel] must be between 1 and 7.");

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
						  if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,
							  MaxPhysCores,true,true,TabLevel[thread_level],i))
						  {
							  poolInterface->DeAllocateAllThreads(true);
							  env->ThrowError("aWarpSharp: Error with the TheadPool while allocating threadpool!");
						  }
						  Offset+=delta;
					  }
				  }
				  else
				  {
					  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					  {
						  poolInterface->DeAllocateAllThreads(true);
						  env->ThrowError("aWarpSharp: Error with the TheadPool while allocating threadpool!");
					  }
				  }
			  }
			  else
			  {
				  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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
	  thread_level=args[10].AsInt(6);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aSobel: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aSobel: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aSobel: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("aSobel: [ThreadLevel] must be between 1 and 7.");

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
						  if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,
							  MaxPhysCores,true,true,TabLevel[thread_level],i))
						  {
							  poolInterface->DeAllocateAllThreads(true);
							  env->ThrowError("aSobel: Error with the TheadPool while allocating threadpool!");
						  }
						  Offset+=delta;
					  }
				  }
				  else
				  {
					  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					  {
						  poolInterface->DeAllocateAllThreads(true);
						  env->ThrowError("aSobel: Error with the TheadPool while allocating threadpool!");
					  }
				  }
			  }
			  else
			  {
				  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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
	  thread_level=args[13].AsInt(6);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aBlur: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aBlur: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aBlur: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("aBlur: [ThreadLevel] must be between 1 and 7.");

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
						  if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,
							  MaxPhysCores,true,true,TabLevel[thread_level],i))
						  {
							  poolInterface->DeAllocateAllThreads(true);
							  env->ThrowError("aBlur: Error with the TheadPool while allocating threadpool!");
						  }
						  Offset+=delta;
					  }
				  }
				  else
				  {
					  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					  {
						  poolInterface->DeAllocateAllThreads(true);
						  env->ThrowError("aBlur: Error with the TheadPool while allocating threadpool!");
					  }
				  }
			  }
			  else
			  {
				  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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
	  thread_level=args[14].AsInt(6);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aWarp: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarp: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarp: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("aWarp: [ThreadLevel] must be between 1 and 7.");

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
						  if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,
							  MaxPhysCores,true,true,TabLevel[thread_level],i))
						  {
							  poolInterface->DeAllocateAllThreads(true);
							  env->ThrowError("aWarp: Error with the TheadPool while allocating threadpool!");
						  }
						  Offset+=delta;
					  }
				  }
				  else
				  {
					  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					  {
						  poolInterface->DeAllocateAllThreads(true);
						  env->ThrowError("aWarp: Error with the TheadPool while allocating threadpool!");
					  }
				  }
			  }
			  else
			  {
				  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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
	  thread_level=args[14].AsInt(6);

	  if (!aWarpSharp_Enable_SSE2) env->ThrowError("aWarp4: SSE2 capable CPU is required");

	  if ((threads<0) || (threads>MAX_MT_THREADS))
		  env->ThrowError("aWarp4: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	  if (prefetch==0) prefetch=1;
	  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		  env->ThrowError("aWarp4: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("aWarp4: [ThreadLevel] must be between 1 and 7.");

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
						  if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,
							  MaxPhysCores,true,true,TabLevel[thread_level],i))
						  {
							  poolInterface->DeAllocateAllThreads(true);
							  env->ThrowError("aWarp4: Error with the TheadPool while allocating threadpool!");
						  }
						  Offset+=delta;
					  }
				  }
				  else
				  {
					  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					  {
						  poolInterface->DeAllocateAllThreads(true);
						  env->ThrowError("aWarp4: Error with the TheadPool while allocating threadpool!");
					  }
				  }
			  }
			  else
			  {
				  if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
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


/***********************/
/*  HDRTools functions */
/***********************/

/*
  Color : int, default value : 2
     0 : BT2100
	 1 : BT2020
	 2 : BT709
	 3 : BT601_525
	 4 : BT601_625
  OutputMode : int, default 0.
     0 : Input 8 Bits -> Output : RGB32, Input > 8 Bits -> Output : RGB64
	 1 : Output : RGB64
	 2 : Output : RGBPS (Planar RGB float)
  HLGMode : bool, default false.
  OOTF : bool, default true.
  EOTF : bool, default true.
  fullrange : bool, default false.
  mpeg2c : bool, default true.
*/
AVSValue __cdecl Create_ConvertYUVtoLinearRGB(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertYUVtoLinearRGB: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if (!vi.IsPlanar() || !vi.IsYUV())
		env->ThrowError("ConvertYUVtoLinearRGB: Input format must be planar YUV");

	const int Color=args[1].AsInt(2);
	int OutputMode=args[2].AsInt(0);
	const uint8_t HDRMode=args[3].AsInt(0);
	const double HLG_Lb=args[4].AsFloat(0.05f);
	const double HLG_Lw=args[5].AsFloat(1000.0f);
	const uint8_t HLGColor=args[6].AsInt(0);
	const bool OOTF=args[7].AsBool(true);
	const bool EOTF=args[8].AsBool(true);
	const bool fullrange=args[9].AsBool(false);
	const bool mpeg2c=args[10].AsBool(true);
	const int threads=args[11].AsInt(0);
	const bool LogicalCores=args[12].AsBool(true);
	const bool MaxPhysCores=args[13].AsBool(true);
	const bool SetAffinity=args[14].AsBool(false);
	const bool sleep = args[15].AsBool(false);
	int prefetch=args[16].AsInt(0);
	int thread_level=args[17].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if (!avsp) OutputMode=0;

	if ((Color<0) || (Color>4))
		env->ThrowError("ConvertYUVtoLinearRGB: [Color] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");
	if ((OutputMode<0) || (OutputMode>2))
		env->ThrowError("ConvertYUVtoLinearRGB: [OutputMode] must be 0, 1 or 2");
	if ((HDRMode<0) || (HDRMode>2))
		env->ThrowError("ConvertYUVtoLinearRGB: [HDRMode] must be 0, 1 or 2");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertYUVtoLinearRGB: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertYUVtoLinearRGB: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertYUVtoLinearRGB: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertYUVtoLinearRGB: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertYUVtoLinearRGB: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertYUVtoLinearRGB: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertYUVtoLinearRGB: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertYUVtoLinearRGB: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertYUVtoLinearRGB(args[0].AsClip(),Color,OutputMode,HDRMode,HLG_Lb,HLG_Lw,HLGColor,OOTF,EOTF,fullrange,
		mpeg2c,threads_number,sleep,env);
}


/*
  Color : int, default value : 2
     0 : BT2100
	 1 : BT2020
	 2 : BT709
	 3 : BT601_525
	 4 : BT601_625
  OutputMode : int, default 0.
     0 : Input 8 Bits -> Output : RGB32, Input > 8 Bits -> Output : RGB64
	 1 : Output : RGB64
	 2 : Output : RGBPS (Planar RGB float)
  HLGMode : bool, default false.
  OOTF : bool, default true.
  EOTF : bool, default true.
  fullrange : bool, default false.
  mpeg2c : bool, default true.
  Rx,Ry,Gx,Gy,Bx,By,Wx,Wy : float, Chromaticity datas.
	Default values are according Color value.
*/
AVSValue __cdecl Create_ConvertYUVtoXYZ(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertYUVtoXYZ: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if (!vi.IsPlanar() || !vi.IsYUV())
		env->ThrowError("ConvertYUVtoXYZ: Input format must be planar YUV");

	float Rx,Ry,Gx,Gy,Bx,By,Wx,Wy;

	const int Color=args[1].AsInt(2);
	int OutputMode=args[2].AsInt(0);
	const uint8_t HDRMode=args[3].AsInt(0);
	const double HLG_Lb=args[4].AsFloat(0.05f);
	const double HLG_Lw=args[5].AsFloat(1000.0f);
	const uint8_t HLGColor=args[6].AsInt(0);
	const bool OOTF=args[7].AsBool(true);
	const bool EOTF=args[8].AsBool(true);
	const bool fullrange=args[9].AsBool(false);
	const bool mpeg2c=args[10].AsBool(true);
	const double Crosstalk=args[19].AsFloat(0.0f);
	const int threads=args[20].AsInt(0);
	const bool LogicalCores=args[21].AsBool(true);
	const bool MaxPhysCores=args[22].AsBool(true);
	const bool SetAffinity=args[23].AsBool(false);
	const bool sleep = args[24].AsBool(false);
	int prefetch=args[25].AsInt(0);
	int thread_level=args[26].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if (!avsp) OutputMode=0;

	if ((Color<0) || (Color>4))
		env->ThrowError("ConvertYUVtoXYZ: [Color] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");
	if ((OutputMode<0) || (OutputMode>2))
		env->ThrowError("ConvertYUVtoXYZ: [OutputMode] must be 0, 1 or 2");
	if ((HDRMode<0) || (HDRMode>2))
		env->ThrowError("ConvertYUVtoXYZ: [HDRMode] must be 0, 1 or 2");
	if ((Crosstalk<0.0) || (Crosstalk>0.33))
		env->ThrowError("ConvertYUVtoXYZ: [Crosstalk] must be in [0..0.33]");

	switch(Color)
	{
		case 0 :
		case 1 :
			Rx=(float)args[11].AsFloat(0.708f);
			Ry=(float)args[12].AsFloat(0.292f);
			Gx=(float)args[13].AsFloat(0.170f);
			Gy=(float)args[14].AsFloat(0.797f);
			Bx=(float)args[15].AsFloat(0.131f);
			By=(float)args[16].AsFloat(0.046f);
			Wx=(float)args[17].AsFloat(0.31271f);
			Wy=(float)args[18].AsFloat(0.32902f);
			break;
		case 2 :
			Rx=(float)args[11].AsFloat(0.640f);
			Ry=(float)args[12].AsFloat(0.330f);
			Gx=(float)args[13].AsFloat(0.300f);
			Gy=(float)args[14].AsFloat(0.600f);
			Bx=(float)args[15].AsFloat(0.150f);
			By=(float)args[16].AsFloat(0.060f);
			Wx=(float)args[17].AsFloat(0.31271f);
			Wy=(float)args[18].AsFloat(0.32902f);
			break;
		case 3 :
			Rx=(float)args[11].AsFloat(0.630f);
			Ry=(float)args[12].AsFloat(0.340f);
			Gx=(float)args[13].AsFloat(0.310f);
			Gy=(float)args[14].AsFloat(0.595f);
			Bx=(float)args[15].AsFloat(0.155f);
			By=(float)args[16].AsFloat(0.070f);
			Wx=(float)args[17].AsFloat(0.31271f);
			Wy=(float)args[18].AsFloat(0.32902f);
			break;
		case 4 :
			Rx=(float)args[11].AsFloat(0.640f);
			Ry=(float)args[12].AsFloat(0.330f);
			Gx=(float)args[13].AsFloat(0.290f);
			Gy=(float)args[14].AsFloat(0.600f);
			Bx=(float)args[15].AsFloat(0.150f);
			By=(float)args[16].AsFloat(0.060f);
			Wx=(float)args[17].AsFloat(0.31271f);
			Wy=(float)args[18].AsFloat(0.32902f);
			break;
	}

	if (((Rx<0.0f) || (Rx>1.0f)) || ((Gx<0.0f) || (Gx>1.0f)) || ((Bx<0.0f) || (Bx>1.0f)) || ((Wx<0.0f) || (Wx>1.0f))
		|| ((Ry<=0.0f) || (Ry>1.0f)) || ((Gy<=0.0f) || (Gy>1.0f)) || ((By<=0.0f) || (By>1.0f)) || ((Wy<=0.0f) || (Wy>1.0f)))
		env->ThrowError("ConvertYUVtoXYZ: Invalid chromaticity datas");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertYUVtoXYZ: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertYUVtoXYZ: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertYUVtoXYZ: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertYUVtoXYZ: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertYUVtoXYZ: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertYUVtoXYZ: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertYUVtoXYZ: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertYUVtoXYZ: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertYUVtoXYZ(args[0].AsClip(),Color,OutputMode,HDRMode,HLG_Lb,HLG_Lw,Crosstalk,HLGColor,OOTF,EOTF,
		fullrange,mpeg2c,Rx,Ry,Gx,Gy,Bx,By,Wx,Wy,threads_number,sleep,env);
}


/*
  Color : int, default value : 2
     0 : BT2100
	 1 : BT2020
	 2 : BT709
	 3 : BT601_525
	 4 : BT601_625
  OutputMode : int, default 0.
     0 : YV24
	 1 : YV16
	 2 : YV12
  HLGMode : bool, default false.
  OOTF : bool, default true.
  OETF : bool, default true.
  fullrange : bool, default false.
  mpeg2c : bool, default true.
  fastmode : bool, default true.
*/
AVSValue __cdecl Create_ConvertLinearRGBtoYUV(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertLinearRGBtoYUV: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR32) && (vi.pixel_type!=VideoInfo::CS_BGR64)
		&& (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertLinearRGBtoYUV: Input format must be RGB32, RGB64 or RGBPS");

	const int Color=args[1].AsInt(2);
	int OutputMode=args[2].AsInt(0);
	const uint8_t HDRMode=args[3].AsInt(0);
	const double HLG_Lb=args[4].AsFloat(0.05f);
	const double HLG_Lw=args[5].AsFloat(1000.0f);
	const uint8_t HLGColor=args[6].AsInt(0);
	const bool OOTF=args[7].AsBool(true);
	const bool EOTF=args[8].AsBool(true);
	const bool fullrange=args[9].AsBool(false);
	const bool mpeg2c=args[10].AsBool(true);
	const bool fastmode=args[11].AsBool(true);
	const int threads=args[12].AsInt(0);
	const bool LogicalCores=args[13].AsBool(true);
	const bool MaxPhysCores=args[14].AsBool(true);
	const bool SetAffinity=args[15].AsBool(false);
	const bool sleep = args[16].AsBool(false);
	int prefetch=args[17].AsInt(0);
	int thread_level=args[18].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((Color<0) || (Color>4))
		env->ThrowError("ConvertLinearRGBtoYUV: [Color] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");
	if ((OutputMode<0) || (OutputMode>2))
		env->ThrowError("ConvertLinearRGBtoYUV: [OutputMode] must be 0, 1 or 2");
	if ((HDRMode<0) || (HDRMode>2))
		env->ThrowError("ConvertLinearRGBtoYUV: [HDRMode] must be 0, 1 or 2");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertLinearRGBtoYUV: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertLinearRGBtoYUV: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertLinearRGBtoYUV: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertLinearRGBtoYUV: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertLinearRGBtoYUV: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertLinearRGBtoYUV: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertLinearRGBtoYUV: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertLinearRGBtoYUV: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertLinearRGBtoYUV(args[0].AsClip(),Color,OutputMode,HDRMode,HLG_Lb,HLG_Lw,HLGColor,OOTF,EOTF,fullrange,
		mpeg2c,fastmode,threads_number,sleep,env);
}


/*
  Color : int, default value : 2
     0 : BT2100
	 1 : BT2020
	 2 : BT709
	 3 : BT601_525
	 4 : BT601_625
  OutputMode : int, default 0.
     0 : No change
	 1 : RGB32 -> RGB64, RGB64 & RGBPS : no change
	 2 : RGB32 & RGB64 -> RGBPS, RGBPS : no change
  HLGMode : bool, default false.
  OOTF : bool, default true.
  EOTF : bool, default true.
  fastmode : bool, default true.
  Rx,Ry,Gx,Gy,Bx,By,Wx,Wy : float, Chromaticity datas.
	Default values are according Color value.
*/
AVSValue __cdecl Create_ConvertRGBtoXYZ(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertRGBtoXYZ: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR32) && (vi.pixel_type!=VideoInfo::CS_BGR64)
		&& (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertRGBtoXYZ: Input format must be Planar RGBPS, RGB32 or RGB64");

	float Rx,Ry,Gx,Gy,Bx,By,Wx,Wy;

	const int Color=args[1].AsInt(2);
	int OutputMode=args[2].AsInt(0);
	const uint8_t HDRMode=args[3].AsInt(0);
	const double HLG_Lb=args[4].AsFloat(0.05f);
	const double HLG_Lw=args[5].AsFloat(1000.0f);
	const uint8_t HLGColor=args[6].AsInt(0);
	const bool OOTF=args[7].AsBool(true);
	const bool EOTF=args[8].AsBool(true);
	const bool fastmode=args[9].AsBool(true);
	const double crosstalk=args[18].AsFloat(0.0f);
	const int threads=args[19].AsInt(0);
	const bool LogicalCores=args[20].AsBool(true);
	const bool MaxPhysCores=args[21].AsBool(true);
	const bool SetAffinity=args[22].AsBool(false);
	const bool sleep = args[23].AsBool(false);
	int prefetch=args[24].AsInt(0);
	int thread_level=args[25].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if (!avsp) OutputMode=0;

	if ((Color<0) || (Color>4))
		env->ThrowError("ConvertRGBtoXYZ: [Color] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");
	if ((OutputMode<0) || (OutputMode>2))
		env->ThrowError("ConvertRGBtoXYZ: [OutputMode] must be 0, 1 or 2");
	if ((HDRMode<0) || (HDRMode>2))
		env->ThrowError("ConvertRGBtoXYZ: [HDRMode] must be 0, 1 or 2");
	if ((crosstalk<0.0) || (crosstalk>0.33))
		env->ThrowError("ConvertRGBtoXYZ: [Crosstalk] must be in [0..0.33]");

	switch(Color)
	{
		case 0 :
		case 1 :
			Rx=(float)args[10].AsFloat(0.708f);
			Ry=(float)args[11].AsFloat(0.292f);
			Gx=(float)args[12].AsFloat(0.170f);
			Gy=(float)args[13].AsFloat(0.797f);
			Bx=(float)args[14].AsFloat(0.131f);
			By=(float)args[15].AsFloat(0.046f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
		case 2 :
			Rx=(float)args[10].AsFloat(0.640f);
			Ry=(float)args[11].AsFloat(0.330f);
			Gx=(float)args[12].AsFloat(0.300f);
			Gy=(float)args[13].AsFloat(0.600f);
			Bx=(float)args[14].AsFloat(0.150f);
			By=(float)args[15].AsFloat(0.060f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
		case 3 :
			Rx=(float)args[10].AsFloat(0.630f);
			Ry=(float)args[11].AsFloat(0.340f);
			Gx=(float)args[12].AsFloat(0.310f);
			Gy=(float)args[13].AsFloat(0.595f);
			Bx=(float)args[14].AsFloat(0.155f);
			By=(float)args[15].AsFloat(0.070f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
		case 4 :
			Rx=(float)args[10].AsFloat(0.640f);
			Ry=(float)args[11].AsFloat(0.330f);
			Gx=(float)args[12].AsFloat(0.290f);
			Gy=(float)args[13].AsFloat(0.600f);
			Bx=(float)args[14].AsFloat(0.150f);
			By=(float)args[15].AsFloat(0.060f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
	}

	if (((Rx<0.0f) || (Rx>1.0f)) || ((Gx<0.0f) || (Gx>1.0f)) || ((Bx<0.0f) || (Bx>1.0f)) || ((Wx<0.0f) || (Wx>1.0f))
		|| ((Ry<=0.0f) || (Ry>1.0f)) || ((Gy<=0.0f) || (Gy>1.0f)) || ((By<=0.0f) || (By>1.0f)) || ((Wy<=0.0f) || (Wy>1.0f)))
		env->ThrowError("ConvertRGBtoXYZ: Invalid chromaticity datas");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertRGBtoXYZ: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertRGBtoXYZ: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertRGBtoXYZ: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertRGBtoXYZ: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertRGBtoXYZ: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertRGBtoXYZ: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertRGBtoXYZ: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertRGBtoXYZ: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertRGBtoXYZ(args[0].AsClip(),Color,OutputMode,HDRMode,HLG_Lb,HLG_Lw,crosstalk,HLGColor,
		OOTF,EOTF,fastmode,Rx,Ry,Gx,Gy,Bx,By,Wx,Wy,threads_number,sleep,env);
}


/*
  Color : int, default value : 2
     0 : BT2100
	 1 : BT2020
	 2 : BT709
	 3 : BT601_525
	 4 : BT601_625
  OutputMode : int, default 0.
     0 : YV24
	 1 : YV16
	 2 : YV12
  HLGMode : bool, default false.
  OOTF : bool, default true.
  OETF : bool, default true.
  fullrange : bool, default false.
  mpeg2c : bool, default true.
  fastmode : bool, default true.
  Rx,Ry,Gx,Gy,Bx,By,Wx,Wy : float, Chromaticity datas.
	Default values are according Color value.
  pColor : int, default value 2 if Color=0, 0 otherwise. Color used in previous YUVtoXYZ.
  pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy : float, Chromaticity datas used in previous YUVtoXYZ.
	Default values are according pColor value.
*/
AVSValue __cdecl Create_ConvertXYZtoYUV(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertXYZtoYUV: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR32) && (vi.pixel_type!=VideoInfo::CS_BGR64)
		&& (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertXYZtoYUV: Input format must be RGB32, RGB64 or RGBPS");

	float Rx,Ry,Gx,Gy,Bx,By,Wx,Wy;
	float pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy;

	const int Color=args[1].AsInt(2);
	int pColor;
	int OutputMode=args[2].AsInt(0);
	const uint8_t HDRMode=args[3].AsInt(0);
	const double HLG_Lb=args[4].AsFloat(0.05f);
	const double HLG_Lw=args[5].AsFloat(1000.0f);
	const uint8_t HLGColor=args[6].AsInt(0);
	const bool OOTF=args[7].AsBool(true);
	const bool EOTF=args[8].AsBool(true);
	const bool fullrange=args[9].AsBool(false);
	const bool mpeg2c=args[10].AsBool(true);
	const bool fastmode=args[11].AsBool(true);
	const double crosstalk=args[29].AsFloat(0.0f);
	const int threads=args[30].AsInt(0);
	const bool LogicalCores=args[31].AsBool(true);
	const bool MaxPhysCores=args[32].AsBool(true);
	const bool SetAffinity=args[33].AsBool(false);
	const bool sleep = args[34].AsBool(false);
	int prefetch=args[35].AsInt(0);
	int thread_level=args[36].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((Color<0) || (Color>4))
		env->ThrowError("ConvertXYZtoYUV: [Color] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");
	switch(Color)
	{
		case 0 : pColor=args[20].AsInt(2); break;
		default : pColor=args[20].AsInt(0); break;
	}
	if ((pColor<0) || (pColor>4))
		env->ThrowError("ConvertXYZtoYUV: [pColor] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");
	if ((OutputMode<0) || (OutputMode>2))
		env->ThrowError("ConvertXYZtoYUV: [OutputMode] must be 0, 1 or 2");
	if ((HDRMode<0) || (HDRMode>2))
		env->ThrowError("ConvertXYZtoYUV: [HDRMode] must be 0, 1 or 2");
	if ((crosstalk<0.0) || (crosstalk>0.33))
		env->ThrowError("ConvertXYZtoYUV: [Crosstalk] must be in [0..0.33]");

	switch(Color)
	{
		case 0 :
		case 1 :
			Rx=(float)args[12].AsFloat(0.708f);
			Ry=(float)args[13].AsFloat(0.292f);
			Gx=(float)args[14].AsFloat(0.170f);
			Gy=(float)args[15].AsFloat(0.797f);
			Bx=(float)args[16].AsFloat(0.131f);
			By=(float)args[17].AsFloat(0.046f);
			Wx=(float)args[18].AsFloat(0.31271f);
			Wy=(float)args[19].AsFloat(0.32902f);
			break;
		case 2 :
			Rx=(float)args[12].AsFloat(0.640f);
			Ry=(float)args[13].AsFloat(0.330f);
			Gx=(float)args[14].AsFloat(0.300f);
			Gy=(float)args[15].AsFloat(0.600f);
			Bx=(float)args[16].AsFloat(0.150f);
			By=(float)args[17].AsFloat(0.060f);
			Wx=(float)args[18].AsFloat(0.31271f);
			Wy=(float)args[19].AsFloat(0.32902f);
			break;
		case 3 :
			Rx=(float)args[12].AsFloat(0.630f);
			Ry=(float)args[13].AsFloat(0.340f);
			Gx=(float)args[14].AsFloat(0.310f);
			Gy=(float)args[15].AsFloat(0.595f);
			Bx=(float)args[16].AsFloat(0.155f);
			By=(float)args[17].AsFloat(0.070f);
			Wx=(float)args[18].AsFloat(0.31271f);
			Wy=(float)args[19].AsFloat(0.32902f);
			break;
		case 4 :
			Rx=(float)args[12].AsFloat(0.640f);
			Ry=(float)args[13].AsFloat(0.330f);
			Gx=(float)args[14].AsFloat(0.290f);
			Gy=(float)args[15].AsFloat(0.600f);
			Bx=(float)args[16].AsFloat(0.150f);
			By=(float)args[17].AsFloat(0.060f);
			Wx=(float)args[18].AsFloat(0.31271f);
			Wy=(float)args[19].AsFloat(0.32902f);
			break;
		default :
			Rx=(float)args[12].AsFloat(0.640f);
			Ry=(float)args[13].AsFloat(0.330f);
			Gx=(float)args[14].AsFloat(0.300f);
			Gy=(float)args[15].AsFloat(0.600f);
			Bx=(float)args[16].AsFloat(0.150f);
			By=(float)args[17].AsFloat(0.060f);
			Wx=(float)args[18].AsFloat(0.31271f);
			Wy=(float)args[19].AsFloat(0.32902f);
			break;
	}

	if (((Rx<0.0f) || (Rx>1.0f)) || ((Gx<0.0f) || (Gx>1.0f)) || ((Bx<0.0f) || (Bx>1.0f)) || ((Wx<0.0f) || (Wx>1.0f))
		|| ((Ry<=0.0f) || (Ry>1.0f)) || ((Gy<=0.0f) || (Gy>1.0f)) || ((By<=0.0f) || (By>1.0f)) || ((Wy<=0.0f) || (Wy>1.0f)))
		env->ThrowError("ConvertXYZtoYUV: Invalid [R,G,B,W][x,y] chromaticity datas");

	switch(pColor)
	{
		case 0 :
		case 1 :
			pRx=(float)args[21].AsFloat(0.708f);
			pRy=(float)args[22].AsFloat(0.292f);
			pGx=(float)args[23].AsFloat(0.170f);
			pGy=(float)args[24].AsFloat(0.797f);
			pBx=(float)args[25].AsFloat(0.131f);
			pBy=(float)args[26].AsFloat(0.046f);
			pWx=(float)args[27].AsFloat(0.31271f);
			pWy=(float)args[28].AsFloat(0.32902f);
			break;
		case 2 :
			pRx=(float)args[21].AsFloat(0.640f);
			pRy=(float)args[22].AsFloat(0.330f);
			pGx=(float)args[23].AsFloat(0.300f);
			pGy=(float)args[24].AsFloat(0.600f);
			pBx=(float)args[25].AsFloat(0.150f);
			pBy=(float)args[26].AsFloat(0.060f);
			pWx=(float)args[27].AsFloat(0.31271f);
			pWy=(float)args[28].AsFloat(0.32902f);
			break;
		case 3 :
			pRx=(float)args[21].AsFloat(0.630f);
			pRy=(float)args[22].AsFloat(0.340f);
			pGx=(float)args[23].AsFloat(0.310f);
			pGy=(float)args[24].AsFloat(0.595f);
			pBx=(float)args[25].AsFloat(0.155f);
			pBy=(float)args[26].AsFloat(0.070f);
			pWx=(float)args[27].AsFloat(0.31271f);
			pWy=(float)args[28].AsFloat(0.32902f);
			break;
		case 4 :
			pRx=(float)args[21].AsFloat(0.640f);
			pRy=(float)args[22].AsFloat(0.330f);
			pGx=(float)args[23].AsFloat(0.290f);
			pGy=(float)args[24].AsFloat(0.600f);
			pBx=(float)args[25].AsFloat(0.150f);
			pBy=(float)args[26].AsFloat(0.060f);
			pWx=(float)args[27].AsFloat(0.31271f);
			pWy=(float)args[28].AsFloat(0.32902f);
			break;
		default :
			pRx=(float)args[21].AsFloat(0.640f);
			pRy=(float)args[22].AsFloat(0.330f);
			pGx=(float)args[23].AsFloat(0.300f);
			pGy=(float)args[24].AsFloat(0.600f);
			pBx=(float)args[25].AsFloat(0.150f);
			pBy=(float)args[26].AsFloat(0.060f);
			pWx=(float)args[27].AsFloat(0.31271f);
			pWy=(float)args[28].AsFloat(0.32902f);
			break;
	}

	if (((pRx<0.0f) || (pRx>1.0f)) || ((pGx<0.0f) || (pGx>1.0f)) || ((pBx<0.0f) || (pBx>1.0f)) || ((pWx<0.0f) || (pWx>1.0f))
		|| ((pRy<=0.0f) || (pRy>1.0f)) || ((pGy<=0.0f) || (pGy>1.0f)) || ((pBy<=0.0f) || (pBy>1.0f)) || ((pWy<=0.0f) || (pWy>1.0f)))
		env->ThrowError("ConvertXYZtoYUV: Invalid [pR,pG,pB,pW][x,y] chromaticity datas");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertXYZtoYUV: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertXYZtoYUV: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertXYZtoYUV: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertXYZtoYUV: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertXYZtoYUV: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertXYZtoYUV: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertXYZtoYUV: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertXYZtoYUV: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertXYZtoYUV(args[0].AsClip(),Color,OutputMode,HDRMode,HLG_Lb,HLG_Lw,crosstalk,HLGColor,
		OOTF,EOTF,fullrange,mpeg2c,fastmode,Rx,Ry,Gx,Gy,Bx,By,Wx,Wy,pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy,
		threads_number,sleep,env);
}


/*
  Color : int, default value : 2
     0 : BT2100
	 1 : BT2020
	 2 : BT709
	 3 : BT601_525
	 4 : BT601_625
  OutputMode : int, default 0.
     0 : No change
	 1 : RGB32 & RGB64 : no change, RGBPS -> RGB64
  HLGMode : bool, default false.
  OOTF : bool, default true.
  OETF : bool, default true.
  fastmode : bool, default true.
  Rx,Ry,Gx,Gy,Bx,By,Wx,Wy : float, Chromaticity datas.
	Default values are according Color value.
  pColor : int, default value 2 if Color=0, 0 otherwise. Color used in previous YUVtoXYZ.
  pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy : float, Chromaticity datas used in previous YUVtoXYZ.
	Default values are according pColor value.
*/
AVSValue __cdecl Create_ConvertXYZtoRGB(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertXYZtoRGB: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR32) && (vi.pixel_type!=VideoInfo::CS_BGR64)
		&& (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertXYZtoRGB: Input format must be RGB32, RGB64 or RGBPS");

	float Rx,Ry,Gx,Gy,Bx,By,Wx,Wy;
	float pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy;

	const int Color=args[1].AsInt(2);
	const int OutputMode=args[2].AsInt(0);
	int pColor;
	const uint8_t HDRMode=args[3].AsInt(0);
	const double HLG_Lb=args[4].AsFloat(0.05f);
	const double HLG_Lw=args[5].AsFloat(1000.0f);
	const uint8_t HLGColor=args[6].AsInt(0);
	const bool OOTF=args[7].AsBool(true);
	const bool EOTF=args[8].AsBool(true);
	const bool fastmode=args[9].AsBool(true);
	const double crosstalk=args[27].AsFloat(0.0f);
	const int threads=args[28].AsInt(0);
	const bool LogicalCores=args[29].AsBool(true);
	const bool MaxPhysCores=args[30].AsBool(true);
	const bool SetAffinity=args[31].AsBool(false);
	const bool sleep = args[32].AsBool(false);
	int prefetch=args[33].AsInt(0);
	int thread_level=args[34].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((Color<0) || (Color>4))
		env->ThrowError("ConvertXYZtoRGB: [Color] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");
	switch(Color)
	{
		case 0 : pColor=args[18].AsInt(2); break;
		default : pColor=args[18].AsInt(0); break;
	}
	if ((pColor<0) || (pColor>4))
		env->ThrowError("ConvertXYZtoRGB: [pColor] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");

	if ((OutputMode<0) || (OutputMode>1))
		env->ThrowError("ConvertXYZtoRGB: [OutputMode] must be 0 or 1");
	if ((HDRMode<0) || (HDRMode>2))
		env->ThrowError("ConvertXYZtoRGB: [HDRMode] must be 0, 1 or 2");
	if ((crosstalk<0.0) || (crosstalk>0.33))
		env->ThrowError("ConvertXYZtoRGB: [Crosstalk] must be in [0..0.33]");

	switch(Color)
	{
		case 0 :
		case 1 :
			Rx=(float)args[10].AsFloat(0.708f);
			Ry=(float)args[11].AsFloat(0.292f);
			Gx=(float)args[12].AsFloat(0.170f);
			Gy=(float)args[13].AsFloat(0.797f);
			Bx=(float)args[14].AsFloat(0.131f);
			By=(float)args[15].AsFloat(0.046f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
		case 2 :
			Rx=(float)args[10].AsFloat(0.640f);
			Ry=(float)args[11].AsFloat(0.330f);
			Gx=(float)args[12].AsFloat(0.300f);
			Gy=(float)args[13].AsFloat(0.600f);
			Bx=(float)args[14].AsFloat(0.150f);
			By=(float)args[15].AsFloat(0.060f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
		case 3 :
			Rx=(float)args[10].AsFloat(0.630f);
			Ry=(float)args[11].AsFloat(0.340f);
			Gx=(float)args[12].AsFloat(0.310f);
			Gy=(float)args[13].AsFloat(0.595f);
			Bx=(float)args[14].AsFloat(0.155f);
			By=(float)args[15].AsFloat(0.070f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
		case 4 :
			Rx=(float)args[10].AsFloat(0.640f);
			Ry=(float)args[11].AsFloat(0.330f);
			Gx=(float)args[12].AsFloat(0.290f);
			Gy=(float)args[13].AsFloat(0.600f);
			Bx=(float)args[14].AsFloat(0.150f);
			By=(float)args[15].AsFloat(0.060f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
		default :
			Rx=(float)args[10].AsFloat(0.640f);
			Ry=(float)args[11].AsFloat(0.330f);
			Gx=(float)args[12].AsFloat(0.300f);
			Gy=(float)args[13].AsFloat(0.600f);
			Bx=(float)args[14].AsFloat(0.150f);
			By=(float)args[15].AsFloat(0.060f);
			Wx=(float)args[16].AsFloat(0.31271f);
			Wy=(float)args[17].AsFloat(0.32902f);
			break;
	}

	if (((Rx<0.0f) || (Rx>1.0f)) || ((Gx<0.0f) || (Gx>1.0f)) || ((Bx<0.0f) || (Bx>1.0f)) || ((Wx<0.0f) || (Wx>1.0f))
		|| ((Ry<=0.0f) || (Ry>1.0f)) || ((Gy<=0.0f) || (Gy>1.0f)) || ((By<=0.0f) || (By>1.0f)) || ((Wy<=0.0f) || (Wy>1.0f)))
		env->ThrowError("ConvertXYZtoRGB: Invalid [R,G,B,W][x,y] chromaticity datas");

	switch(pColor)
	{
		case 0 :
		case 1 :
			pRx=(float)args[19].AsFloat(0.708f);
			pRy=(float)args[20].AsFloat(0.292f);
			pGx=(float)args[21].AsFloat(0.170f);
			pGy=(float)args[22].AsFloat(0.797f);
			pBx=(float)args[23].AsFloat(0.131f);
			pBy=(float)args[24].AsFloat(0.046f);
			pWx=(float)args[25].AsFloat(0.31271f);
			pWy=(float)args[26].AsFloat(0.32902f);
			break;
		case 2 :
			pRx=(float)args[19].AsFloat(0.640f);
			pRy=(float)args[20].AsFloat(0.330f);
			pGx=(float)args[21].AsFloat(0.300f);
			pGy=(float)args[22].AsFloat(0.600f);
			pBx=(float)args[23].AsFloat(0.150f);
			pBy=(float)args[24].AsFloat(0.060f);
			pWx=(float)args[25].AsFloat(0.31271f);
			pWy=(float)args[26].AsFloat(0.32902f);
			break;
		case 3 :
			pRx=(float)args[19].AsFloat(0.630f);
			pRy=(float)args[20].AsFloat(0.340f);
			pGx=(float)args[21].AsFloat(0.310f);
			pGy=(float)args[22].AsFloat(0.595f);
			pBx=(float)args[23].AsFloat(0.155f);
			pBy=(float)args[24].AsFloat(0.070f);
			pWx=(float)args[25].AsFloat(0.31271f);
			pWy=(float)args[26].AsFloat(0.32902f);
			break;
		case 4 :
			pRx=(float)args[19].AsFloat(0.640f);
			pRy=(float)args[20].AsFloat(0.330f);
			pGx=(float)args[21].AsFloat(0.290f);
			pGy=(float)args[22].AsFloat(0.600f);
			pBx=(float)args[23].AsFloat(0.150f);
			pBy=(float)args[24].AsFloat(0.060f);
			pWx=(float)args[25].AsFloat(0.31271f);
			pWy=(float)args[26].AsFloat(0.32902f);
			break;
		default :
			pRx=(float)args[19].AsFloat(0.640f);
			pRy=(float)args[20].AsFloat(0.330f);
			pGx=(float)args[21].AsFloat(0.300f);
			pGy=(float)args[22].AsFloat(0.600f);
			pBx=(float)args[23].AsFloat(0.150f);
			pBy=(float)args[24].AsFloat(0.060f);
			pWx=(float)args[25].AsFloat(0.31271f);
			pWy=(float)args[26].AsFloat(0.32902f);
			break;
	}

	if (((pRx<0.0f) || (pRx>1.0f)) || ((pGx<0.0f) || (pGx>1.0f)) || ((pBx<0.0f) || (pBx>1.0f)) || ((pWx<0.0f) || (pWx>1.0f))
		|| ((pRy<=0.0f) || (pRy>1.0f)) || ((pGy<=0.0f) || (pGy>1.0f)) || ((pBy<=0.0f) || (pBy>1.0f)) || ((pWy<=0.0f) || (pWy>1.0f)))
		env->ThrowError("ConvertXYZtoRGB: Invalid [pR,pG,pB,pW][x,y] chromaticity datas");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertXYZtoRGB: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertXYZtoRGB: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertXYZtoRGB: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertXYZtoRGB: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertXYZtoRGB: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertXYZtoRGB: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertXYZtoRGB: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertXYZtoRGB: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertXYZtoRGB(args[0].AsClip(),Color,OutputMode,HDRMode,HLG_Lb,HLG_Lw,crosstalk,HLGColor,
		OOTF,EOTF,fastmode,Rx,Ry,Gx,Gy,Bx,By,Wx,Wy,pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy,
		threads_number,sleep,env);
}


/*
  MinMastering,MaxMastering : SEI data if avaible from video. Default 0.0,1000.0.
*/
AVSValue __cdecl Create_ConvertXYZ_Scale_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: Input format must be RGB64 or RGBPS");

	const float Coeff_X=(float)args[1].AsFloat(100.0f);
	const float Coeff_Y=(float)args[2].AsFloat(Coeff_X);
	const float Coeff_Z=(float)args[3].AsFloat(Coeff_X);
	const int threads=args[4].AsInt(0);
	const bool LogicalCores=args[5].AsBool(true);
	const bool MaxPhysCores=args[6].AsBool(true);
	const bool SetAffinity=args[7].AsBool(false);
	const bool sleep = args[8].AsBool(false);
	int prefetch=args[9].AsInt(0);
	int thread_level=args[10].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertXYZ_Scale_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertXYZ_Scale_HDRtoSDR(args[0].AsClip(),Coeff_X,Coeff_Y,Coeff_Z,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertXYZ_Scale_SDRtoHDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: Input format must be RGB64 or RGBPS");

	const float Coeff_X=(float)args[1].AsFloat(100.0f);
	const float Coeff_Y=(float)args[2].AsFloat(Coeff_X);
	const float Coeff_Z=(float)args[3].AsFloat(Coeff_X);
	const int threads=args[4].AsInt(0);
	const bool LogicalCores=args[5].AsBool(true);
	const bool MaxPhysCores=args[6].AsBool(true);
	const bool SetAffinity=args[7].AsBool(false);
	const bool sleep = args[8].AsBool(false);
	int prefetch=args[9].AsInt(0);
	int thread_level=args[10].AsInt(6);

	if ((Coeff_X==0.0f) || (Coeff_Y==0.0f) || (Coeff_Z==0.0f))
		env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: Wrong parameter value!");

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertXYZ_Scale_SDRtoHDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertXYZ_Scale_SDRtoHDR(args[0].AsClip(),Coeff_X,Coeff_Y,Coeff_Z,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertXYZ_Hable_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Input format must be RGB64 or RGBPS");


	const double exp_X=args[1].AsFloat(2.0f);
	const double w_X=args[2].AsFloat(11.2f);
	const double a_X=args[3].AsFloat(0.15f);
	const double b_X=args[4].AsFloat(0.5f);
	const double c_X=args[5].AsFloat(0.1f);
	const double d_X=args[6].AsFloat(0.2f);
	const double e_X=args[7].AsFloat(0.02f);
	const double f_X=args[8].AsFloat(0.3f);

	const double exp_Y=args[9].AsFloat((float)exp_X);
	const double w_Y=args[10].AsFloat((float)w_X);
	const double a_Y=args[11].AsFloat((float)a_X);
	const double b_Y=args[12].AsFloat((float)b_X);
	const double c_Y=args[13].AsFloat((float)c_X);
	const double d_Y=args[14].AsFloat((float)d_X);
	const double e_Y=args[15].AsFloat((float)e_X);
	const double f_Y=args[16].AsFloat((float)f_X);

	const double exp_Z=args[17].AsFloat((float)exp_X);
	const double w_Z=args[18].AsFloat((float)w_X);
	const double a_Z=args[19].AsFloat((float)a_X);
	const double b_Z=args[20].AsFloat((float)b_X);
	const double c_Z=args[21].AsFloat((float)c_X);
	const double d_Z=args[22].AsFloat((float)d_X);
	const double e_Z=args[23].AsFloat((float)e_X);
	const double f_Z=args[24].AsFloat((float)f_X);

	if ((f_X==0.0) || (f_Y==0.0) || (f_Z==0.0))
		env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Wrong parameter value!");

	const uint8_t pColor=args[25].AsInt(0);

	if ((pColor<0) || (pColor>4))
		env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: [pColor] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");

	float pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy;

	switch(pColor)
	{
		case 0 :
		case 1 :
			pRx=(float)args[26].AsFloat(0.708f);
			pRy=(float)args[27].AsFloat(0.292f);
			pGx=(float)args[28].AsFloat(0.170f);
			pGy=(float)args[29].AsFloat(0.797f);
			pBx=(float)args[30].AsFloat(0.131f);
			pBy=(float)args[31].AsFloat(0.046f);
			pWx=(float)args[32].AsFloat(0.31271f);
			pWy=(float)args[33].AsFloat(0.32902f);
			break;
		case 2 :
			pRx=(float)args[26].AsFloat(0.640f);
			pRy=(float)args[27].AsFloat(0.330f);
			pGx=(float)args[28].AsFloat(0.300f);
			pGy=(float)args[29].AsFloat(0.600f);
			pBx=(float)args[30].AsFloat(0.150f);
			pBy=(float)args[31].AsFloat(0.060f);
			pWx=(float)args[32].AsFloat(0.31271f);
			pWy=(float)args[33].AsFloat(0.32902f);
			break;
		case 3 :
			pRx=(float)args[26].AsFloat(0.630f);
			pRy=(float)args[27].AsFloat(0.340f);
			pGx=(float)args[28].AsFloat(0.310f);
			pGy=(float)args[29].AsFloat(0.595f);
			pBx=(float)args[30].AsFloat(0.155f);
			pBy=(float)args[31].AsFloat(0.070f);
			pWx=(float)args[32].AsFloat(0.31271f);
			pWy=(float)args[33].AsFloat(0.32902f);
			break;
		case 4 :
			pRx=(float)args[26].AsFloat(0.640f);
			pRy=(float)args[27].AsFloat(0.330f);
			pGx=(float)args[28].AsFloat(0.290f);
			pGy=(float)args[29].AsFloat(0.600f);
			pBx=(float)args[30].AsFloat(0.150f);
			pBy=(float)args[31].AsFloat(0.060f);
			pWx=(float)args[32].AsFloat(0.31271f);
			pWy=(float)args[33].AsFloat(0.32902f);
			break;
		default :
			pRx=(float)args[26].AsFloat(0.640f);
			pRy=(float)args[27].AsFloat(0.330f);
			pGx=(float)args[28].AsFloat(0.300f);
			pGy=(float)args[29].AsFloat(0.600f);
			pBx=(float)args[30].AsFloat(0.150f);
			pBy=(float)args[31].AsFloat(0.060f);
			pWx=(float)args[32].AsFloat(0.31271f);
			pWy=(float)args[33].AsFloat(0.32902f);
			break;
	}

	if (((pRx<0.0f) || (pRx>1.0f)) || ((pGx<0.0f) || (pGx>1.0f)) || ((pBx<0.0f) || (pBx>1.0f)) || ((pWx<0.0f) || (pWx>1.0f))
		|| ((pRy<=0.0f) || (pRy>1.0f)) || ((pGy<=0.0f) || (pGy>1.0f)) || ((pBy<=0.0f) || (pBy>1.0f)) || ((pWy<=0.0f) || (pWy>1.0f)))
		env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Invalid [pR,pG,pB,pW][x,y] chromaticity datas");

	const bool fastmode=args[34].AsBool(true);

	const int threads=args[35].AsInt(0);
	const bool LogicalCores=args[36].AsBool(true);
	const bool MaxPhysCores=args[37].AsBool(true);
	const bool SetAffinity=args[38].AsBool(false);
	const bool sleep = args[39].AsBool(false);
	int prefetch=args[40].AsInt(0);
	int thread_level=args[41].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertXYZ_Hable_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertXYZ_Hable_HDRtoSDR(args[0].AsClip(),exp_X,w_X,a_X,b_X,c_X,d_X,e_X,f_X,exp_Y,w_Y,a_Y,b_Y,c_Y,d_Y,e_Y,f_Y,
		exp_Z,w_Z,a_Z,b_Z,c_Z,d_Z,e_Z,f_Z,pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy,fastmode,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertRGB_Hable_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertRGB_Hable_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertRGB_Hable_HDRtoSDR: Input format must be RGB64 or RGBPS");


	const double exp_R=args[1].AsFloat(2.0f);
	const double w_R=args[2].AsFloat(11.2f);
	const double a_R=args[3].AsFloat(0.15f);
	const double b_R=args[4].AsFloat(0.5f);
	const double c_R=args[5].AsFloat(0.1f);
	const double d_R=args[6].AsFloat(0.2f);
	const double e_R=args[7].AsFloat(0.02f);
	const double f_R=args[8].AsFloat(0.3f);

	const double exp_G=args[9].AsFloat((float)exp_R);
	const double w_G=args[10].AsFloat((float)w_R);
	const double a_G=args[11].AsFloat((float)a_R);
	const double b_G=args[12].AsFloat((float)b_R);
	const double c_G=args[13].AsFloat((float)c_R);
	const double d_G=args[14].AsFloat((float)d_R);
	const double e_G=args[15].AsFloat((float)e_R);
	const double f_G=args[16].AsFloat((float)f_R);

	const double exp_B=args[17].AsFloat((float)exp_R);
	const double w_B=args[18].AsFloat((float)w_R);
	const double a_B=args[19].AsFloat((float)a_R);
	const double b_B=args[20].AsFloat((float)b_R);
	const double c_B=args[21].AsFloat((float)c_R);
	const double d_B=args[22].AsFloat((float)d_R);
	const double e_B=args[23].AsFloat((float)e_R);
	const double f_B=args[24].AsFloat((float)f_R);

	if ((f_R==0.0) || (f_G==0.0) || (f_B==0.0))
		env->ThrowError("ConvertRGB_Hable_HDRtoSDR: Wrong parameter value!");

	const bool fastmode=args[25].AsBool(true);

	const int threads=args[26].AsInt(0);
	const bool LogicalCores=args[27].AsBool(true);
	const bool MaxPhysCores=args[28].AsBool(true);
	const bool SetAffinity=args[29].AsBool(false);
	const bool sleep = args[30].AsBool(false);
	int prefetch=args[31].AsInt(0);
	int thread_level=args[32].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertRGB_Hable_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertRGB_Hable_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertRGB_Hable_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertRGB_Hable_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertRGB_Hable_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertRGB_Hable_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertRGB_Hable_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertRGB_Hable_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertRGB_Hable_HDRtoSDR(args[0].AsClip(),exp_R,w_R,a_R,b_R,c_R,d_R,e_R,f_R,exp_G,w_G,a_G,b_G,c_G,d_G,e_G,f_G,
		exp_B,w_B,a_B,b_B,c_B,d_B,e_B,f_B,fastmode,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertXYZ_Mobius_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: Input format must be RGB64 or RGBPS");


	const double exp_X=args[1].AsFloat(2.0f);
	const double trans_X=args[2].AsFloat(0.3f);
	const double peak_X=args[3].AsFloat(1.0f);

	const double exp_Y=args[4].AsFloat((float)exp_X);
	const double trans_Y=args[5].AsFloat((float)trans_X);
	const double peak_Y=args[6].AsFloat((float)peak_X);

	const double exp_Z=args[7].AsFloat((float)exp_X);
	const double trans_Z=args[8].AsFloat((float)trans_X);
	const double peak_Z=args[9].AsFloat((float)peak_X);

	const uint8_t pColor=args[10].AsInt(0);

	if ((pColor<0) || (pColor>4))
		env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: [pColor] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");

	float pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy;

	switch(pColor)
	{
		case 0 :
		case 1 :
			pRx=(float)args[11].AsFloat(0.708f);
			pRy=(float)args[12].AsFloat(0.292f);
			pGx=(float)args[13].AsFloat(0.170f);
			pGy=(float)args[14].AsFloat(0.797f);
			pBx=(float)args[15].AsFloat(0.131f);
			pBy=(float)args[16].AsFloat(0.046f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		case 2 :
			pRx=(float)args[11].AsFloat(0.640f);
			pRy=(float)args[12].AsFloat(0.330f);
			pGx=(float)args[13].AsFloat(0.300f);
			pGy=(float)args[14].AsFloat(0.600f);
			pBx=(float)args[15].AsFloat(0.150f);
			pBy=(float)args[16].AsFloat(0.060f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		case 3 :
			pRx=(float)args[11].AsFloat(0.630f);
			pRy=(float)args[12].AsFloat(0.340f);
			pGx=(float)args[13].AsFloat(0.310f);
			pGy=(float)args[14].AsFloat(0.595f);
			pBx=(float)args[15].AsFloat(0.155f);
			pBy=(float)args[16].AsFloat(0.070f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		case 4 :
			pRx=(float)args[11].AsFloat(0.640f);
			pRy=(float)args[12].AsFloat(0.330f);
			pGx=(float)args[13].AsFloat(0.290f);
			pGy=(float)args[14].AsFloat(0.600f);
			pBx=(float)args[15].AsFloat(0.150f);
			pBy=(float)args[16].AsFloat(0.060f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		default :
			pRx=(float)args[11].AsFloat(0.640f);
			pRy=(float)args[12].AsFloat(0.330f);
			pGx=(float)args[13].AsFloat(0.300f);
			pGy=(float)args[14].AsFloat(0.600f);
			pBx=(float)args[15].AsFloat(0.150f);
			pBy=(float)args[16].AsFloat(0.060f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
	}

	if (((pRx<0.0f) || (pRx>1.0f)) || ((pGx<0.0f) || (pGx>1.0f)) || ((pBx<0.0f) || (pBx>1.0f)) || ((pWx<0.0f) || (pWx>1.0f))
		|| ((pRy<=0.0f) || (pRy>1.0f)) || ((pGy<=0.0f) || (pGy>1.0f)) || ((pBy<=0.0f) || (pBy>1.0f)) || ((pWy<=0.0f) || (pWy>1.0f)))
		env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: Invalid [pR,pG,pB,pW][x,y] chromaticity datas");

	const bool fastmode=args[19].AsBool(true);

	const int threads=args[20].AsInt(0);
	const bool LogicalCores=args[21].AsBool(true);
	const bool MaxPhysCores=args[22].AsBool(true);
	const bool SetAffinity=args[23].AsBool(false);
	const bool sleep = args[24].AsBool(false);
	int prefetch=args[25].AsInt(0);
	int thread_level=args[26].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertXYZ_Mobius_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertXYZ_Mobius_HDRtoSDR(args[0].AsClip(),exp_X,trans_X,peak_X,exp_Y,trans_Y,peak_Y,
		exp_Z,trans_Z,peak_Z,pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy,fastmode,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertRGB_Mobius_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: Input format must be RGB64 or RGBPS");


	const double exp_R=args[1].AsFloat(2.0f);
	const double trans_R=args[2].AsFloat(0.3f);
	const double peak_R=args[3].AsFloat(1.0f);

	const double exp_G=args[4].AsFloat((float)exp_R);
	const double trans_G=args[5].AsFloat((float)trans_R);
	const double peak_G=args[6].AsFloat((float)peak_R);

	const double exp_B=args[7].AsFloat((float)exp_R);
	const double trans_B=args[8].AsFloat((float)trans_R);
	const double peak_B=args[9].AsFloat((float)peak_R);

	const bool fastmode=args[10].AsBool(true);

	const int threads=args[11].AsInt(0);
	const bool LogicalCores=args[12].AsBool(true);
	const bool MaxPhysCores=args[13].AsBool(true);
	const bool SetAffinity=args[14].AsBool(false);
	const bool sleep = args[15].AsBool(false);
	int prefetch=args[16].AsInt(0);
	int thread_level=args[17].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertRGB_Mobius_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertRGB_Mobius_HDRtoSDR(args[0].AsClip(),exp_R,trans_R,peak_R,exp_G,trans_G,peak_G,
		exp_B,trans_B,peak_B,fastmode,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertXYZ_Reinhard_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Input format must be RGB64 or RGBPS");


	const double exp_X=args[1].AsFloat(1.5f);
	const double contr_X=args[2].AsFloat(0.5f);
	const double peak_X=args[3].AsFloat(1.0f);

	const double exp_Y=args[4].AsFloat((float)exp_X);
	const double contr_Y=args[5].AsFloat((float)contr_X);
	const double peak_Y=args[6].AsFloat((float)peak_X);

	const double exp_Z=args[7].AsFloat((float)exp_X);
	const double contr_Z=args[8].AsFloat((float)contr_X);
	const double peak_Z=args[9].AsFloat((float)peak_X);
		
	if ((contr_X==0.0) || (contr_Y==0.0) || (contr_Z==0.0) ||
		(peak_X==0.0) || (peak_Y==0.0) || (peak_Z==0.0))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Wrong parameter value!");

	const uint8_t pColor=args[10].AsInt(0);

	if ((pColor<0) || (pColor>4))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: [pColor] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");

	float pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy;

	switch(pColor)
	{
		case 0 :
		case 1 :
			pRx=(float)args[11].AsFloat(0.708f);
			pRy=(float)args[12].AsFloat(0.292f);
			pGx=(float)args[13].AsFloat(0.170f);
			pGy=(float)args[14].AsFloat(0.797f);
			pBx=(float)args[15].AsFloat(0.131f);
			pBy=(float)args[16].AsFloat(0.046f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		case 2 :
			pRx=(float)args[11].AsFloat(0.640f);
			pRy=(float)args[12].AsFloat(0.330f);
			pGx=(float)args[13].AsFloat(0.300f);
			pGy=(float)args[14].AsFloat(0.600f);
			pBx=(float)args[15].AsFloat(0.150f);
			pBy=(float)args[16].AsFloat(0.060f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		case 3 :
			pRx=(float)args[11].AsFloat(0.630f);
			pRy=(float)args[12].AsFloat(0.340f);
			pGx=(float)args[13].AsFloat(0.310f);
			pGy=(float)args[14].AsFloat(0.595f);
			pBx=(float)args[15].AsFloat(0.155f);
			pBy=(float)args[16].AsFloat(0.070f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		case 4 :
			pRx=(float)args[11].AsFloat(0.640f);
			pRy=(float)args[12].AsFloat(0.330f);
			pGx=(float)args[13].AsFloat(0.290f);
			pGy=(float)args[14].AsFloat(0.600f);
			pBx=(float)args[15].AsFloat(0.150f);
			pBy=(float)args[16].AsFloat(0.060f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
		default :
			pRx=(float)args[11].AsFloat(0.640f);
			pRy=(float)args[12].AsFloat(0.330f);
			pGx=(float)args[13].AsFloat(0.300f);
			pGy=(float)args[14].AsFloat(0.600f);
			pBx=(float)args[15].AsFloat(0.150f);
			pBy=(float)args[16].AsFloat(0.060f);
			pWx=(float)args[17].AsFloat(0.31271f);
			pWy=(float)args[18].AsFloat(0.32902f);
			break;
	}

	if (((pRx<0.0f) || (pRx>1.0f)) || ((pGx<0.0f) || (pGx>1.0f)) || ((pBx<0.0f) || (pBx>1.0f)) || ((pWx<0.0f) || (pWx>1.0f))
		|| ((pRy<=0.0f) || (pRy>1.0f)) || ((pGy<=0.0f) || (pGy>1.0f)) || ((pBy<=0.0f) || (pBy>1.0f)) || ((pWy<=0.0f) || (pWy>1.0f)))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Invalid [pR,pG,pB,pW][x,y] chromaticity datas");

	const bool fastmode=args[19].AsBool(true);

	const int threads=args[20].AsInt(0);
	const bool LogicalCores=args[21].AsBool(true);
	const bool MaxPhysCores=args[22].AsBool(true);
	const bool SetAffinity=args[23].AsBool(false);
	const bool sleep = args[24].AsBool(false);
	int prefetch=args[25].AsInt(0);
	int thread_level=args[26].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertXYZ_Reinhard_HDRtoSDR(args[0].AsClip(),exp_X,contr_X,peak_X,exp_Y,contr_Y,peak_Y,
		exp_Z,contr_Z,peak_Z,pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy,fastmode,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertRGB_Reinhard_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: Input format must be RGB64 or RGBPS");


	const double exp_R=args[1].AsFloat(1.5f);
	const double contr_R=args[2].AsFloat(0.5f);
	const double peak_R=args[3].AsFloat(1.0f);

	const double exp_G=args[4].AsFloat((float)exp_R);
	const double contr_G=args[5].AsFloat((float)contr_R);
	const double peak_G=args[6].AsFloat((float)peak_R);

	const double exp_B=args[7].AsFloat((float)exp_R);
	const double contr_B=args[8].AsFloat((float)contr_R);
	const double peak_B=args[9].AsFloat((float)peak_R);

	const bool fastmode=args[10].AsBool(true);

	const int threads=args[11].AsInt(0);
	const bool LogicalCores=args[12].AsBool(true);
	const bool MaxPhysCores=args[13].AsBool(true);
	const bool SetAffinity=args[14].AsBool(false);
	const bool sleep = args[15].AsBool(false);
	int prefetch=args[16].AsInt(0);
	int thread_level=args[17].AsInt(6);

	if ((contr_R==0.0) || (contr_G==0.0) || (contr_B==0.0) ||
		(peak_R==0.0) || (peak_G==0.0) || (peak_B==0.0))
		env->ThrowError("ConvertXYZ_Reinhard_HDRtoSDR: Wrong parameter value!");

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertRGB_Reinhard_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertRGB_Reinhard_HDRtoSDR(args[0].AsClip(),exp_R,contr_R,peak_R,exp_G,contr_G,peak_G,
		exp_B,contr_B,peak_B,fastmode,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: Input format must be RGB64 or RGBPS");

	const double Lhdr=args[1].AsFloat(1000.0f);
	const double Lsdr=args[2].AsFloat(100.0f);
	const double CoeffAdj=args[3].AsFloat(1.0f);

	const bool fastmode=args[4].AsBool(true);

	const int threads=args[5].AsInt(0);
	const bool LogicalCores=args[6].AsBool(true);
	const bool MaxPhysCores=args[7].AsBool(true);
	const bool SetAffinity=args[8].AsBool(false);
	const bool sleep=args[9].AsBool(false);
	int prefetch=args[10].AsInt(0);
	int thread_level=args[11].AsInt(6);

	if ((Lhdr<=0.0) || (Lhdr>10000.0) || (Lsdr<=0.0) || (Lsdr>100.0))
		env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: Wrong parameter value!");

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR(args[0].AsClip(),Lhdr,Lsdr,CoeffAdj,
		fastmode,threads_number,sleep,env);
}


AVSValue __cdecl Create_ConverXYZ_BT2446_C_HDRtoSDR(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip()) env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: arg 0 must be a clip !");

	VideoInfo vi = args[0].AsClip()->GetVideoInfo();

	if ((vi.pixel_type!=VideoInfo::CS_BGR64) && (vi.pixel_type!=VideoInfo::CS_RGBPS))
		env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: Input format must be RGB64 or RGBPS");

	const bool ChromaC=args[1].AsBool(false);
	const bool PQMode=args[2].AsBool(false);
	const float Lhdr=(PQMode)?(float)args[3].AsFloat(10000.0f):(float)args[3].AsFloat(1000.0f);
	const float Lsdr=(float)args[4].AsFloat(100.0f);
	const float pct_ref=(PQMode)?(float)args[5].AsFloat(0.58f):(float)args[5].AsFloat(0.75f);
	const float pct_ip=(float)args[6].AsFloat(0.80f);
	const float pct_wp=(float)args[7].AsFloat(0.96f);
	const float pct_sdr_skin=(float)args[8].AsFloat(0.70f);
	const float pct_hdr_skin=(PQMode)?(float)args[9].AsFloat(0.44f):(float)args[9].AsFloat(0.50f);
	const float WhiteShift=(float)args[10].AsFloat(0.0f);

	if (((pct_ref<0) || (pct_ref>1.0)) || ((pct_ip<0) || (pct_ip>1.0)) || ((pct_wp<0) || (pct_wp>1.0))
		|| ((pct_sdr_skin<0) || (pct_sdr_skin>1.0)) || ((pct_hdr_skin<0) || (pct_hdr_skin>1.0)))
		env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: [pct_xx] must be in [0.0..1.0]");

	const uint8_t pColor=args[11].AsInt(0);

	if ((pColor<0) || (pColor>4))
		env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: [pColor] must be 0 (BT2100), 1 (BT2020), 2 (BT709), 3 (BT601_525), 4 (BT601_625)");

	float pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy;

	switch(pColor)
	{
		case 0 :
		case 1 :
			pRx=(float)args[12].AsFloat(0.708f);
			pRy=(float)args[13].AsFloat(0.292f);
			pGx=(float)args[14].AsFloat(0.170f);
			pGy=(float)args[15].AsFloat(0.797f);
			pBx=(float)args[16].AsFloat(0.131f);
			pBy=(float)args[17].AsFloat(0.046f);
			pWx=(float)args[18].AsFloat(0.31271f);
			pWy=(float)args[19].AsFloat(0.32902f);
			break;
		case 2 :
			pRx=(float)args[12].AsFloat(0.640f);
			pRy=(float)args[13].AsFloat(0.330f);
			pGx=(float)args[14].AsFloat(0.300f);
			pGy=(float)args[15].AsFloat(0.600f);
			pBx=(float)args[16].AsFloat(0.150f);
			pBy=(float)args[17].AsFloat(0.060f);
			pWx=(float)args[18].AsFloat(0.31271f);
			pWy=(float)args[19].AsFloat(0.32902f);
			break;
		case 3 :
			pRx=(float)args[12].AsFloat(0.630f);
			pRy=(float)args[13].AsFloat(0.340f);
			pGx=(float)args[14].AsFloat(0.310f);
			pGy=(float)args[15].AsFloat(0.595f);
			pBx=(float)args[16].AsFloat(0.155f);
			pBy=(float)args[17].AsFloat(0.070f);
			pWx=(float)args[18].AsFloat(0.31271f);
			pWy=(float)args[19].AsFloat(0.32902f);
			break;
		case 4 :
			pRx=(float)args[12].AsFloat(0.640f);
			pRy=(float)args[13].AsFloat(0.330f);
			pGx=(float)args[14].AsFloat(0.290f);
			pGy=(float)args[15].AsFloat(0.600f);
			pBx=(float)args[16].AsFloat(0.150f);
			pBy=(float)args[17].AsFloat(0.060f);
			pWx=(float)args[18].AsFloat(0.31271f);
			pWy=(float)args[19].AsFloat(0.32902f);
			break;
		default :
			pRx=(float)args[12].AsFloat(0.640f);
			pRy=(float)args[13].AsFloat(0.330f);
			pGx=(float)args[14].AsFloat(0.300f);
			pGy=(float)args[15].AsFloat(0.600f);
			pBx=(float)args[16].AsFloat(0.150f);
			pBy=(float)args[17].AsFloat(0.060f);
			pWx=(float)args[18].AsFloat(0.31271f);
			pWy=(float)args[19].AsFloat(0.32902f);
			break;
	}

	if (((pRx<0.0f) || (pRx>1.0f)) || ((pGx<0.0f) || (pGx>1.0f)) || ((pBx<0.0f) || (pBx>1.0f)) || ((pWx<0.0f) || (pWx>1.0f))
		|| ((pRy<=0.0f) || (pRy>1.0f)) || ((pGy<=0.0f) || (pGy>1.0f)) || ((pBy<=0.0f) || (pBy>1.0f)) || ((pWy<=0.0f) || (pWy>1.0f)))
		env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: Invalid [pR,pG,pB,pW][x,y] chromaticity datas");


	const bool fastmode=args[20].AsBool(true);

	const int threads=args[21].AsInt(0);
	const bool LogicalCores=args[22].AsBool(true);
	const bool MaxPhysCores=args[23].AsBool(true);
	const bool SetAffinity=args[24].AsBool(false);
	const bool sleep=args[25].AsBool(false);
	int prefetch=args[26].AsInt(0);
	int thread_level=args[27].AsInt(6);

	const bool avsp=env->FunctionExists("ConvertBits");

	if ((threads<0) || (threads>MAX_MT_THREADS))
		env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: [threads] must be between 0 and %ld.",MAX_MT_THREADS);
	if (prefetch==0) prefetch=1;
	if ((prefetch<0) || (prefetch>MAX_THREAD_POOL))
		env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
	if ((thread_level<1) || (thread_level>7))
		env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: [ThreadLevel] must be between 1 and 7.");

	uint8_t threads_number=1;

	if (threads!=1)
	{
		const ThreadLevelName TabLevel[8]={NoneThreadLevel,IdleThreadLevel,LowestThreadLevel,
			BelowThreadLevel,NormalThreadLevel,AboveThreadLevel,HighestThreadLevel,CriticalThreadLevel};

		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(threads,LogicalCores);

		if (threads_number==0) env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (prefetch>1)
			{
				if (SetAffinity && (prefetch<=poolInterface->GetPhysicalCoreNumber()))
				{
					float delta=(float)poolInterface->GetPhysicalCoreNumber()/(float)prefetch,Offset=0.0f;

					for(uint8_t i=0; i<prefetch; i++)
					{
						if (!poolInterface->AllocateThreads(threads_number,(uint8_t)ceil(Offset),0,MaxPhysCores,
							true,true,TabLevel[thread_level],i))
						{
							poolInterface->DeAllocateAllThreads(true);
							env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
						}
						Offset+=delta;
					}
				}
				else
				{
					if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,false,true,TabLevel[thread_level],-1))
					{
						poolInterface->DeAllocateAllThreads(true);
						env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
					}
				}
			}
			else
			{
				if (!poolInterface->AllocateThreads(threads_number,0,0,MaxPhysCores,SetAffinity,true,TabLevel[thread_level],-1))
				{
					poolInterface->DeAllocateAllThreads(true);
					env->ThrowError("ConverXYZ_BT2446_C_HDRtoSDR: Error with the TheadPool while allocating threadpool!");
				}
			}
		}
	}

	return new ConverXYZ_BT2446_C_HDRtoSDR(args[0].AsClip(),ChromaC,PQMode,Lhdr,Lsdr,pct_ref,pct_ip,pct_wp,
		pct_sdr_skin,pct_hdr_skin,WhiteShift,pRx,pRy,pGx,pGy,pBx,pBy,pWx,pWy,fastmode,threads_number,sleep,env);
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
		"c[threshold]i[mode]i[output]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_AutoYUY2, 0);

	env->AddFunction("nnedi3", "c[field]i[dh]b[Y]b[U]b[V]b[nsize]i[nns]i[qual]i[etype]i[pscrn]i" \
		"[threads]i[opt]i[fapprox]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[A]b[sleep]b[prefetch]i" \
		"[range]i[ThreadLevel]i", Create_nnedi3, 0);
	env->AddFunction("nnedi3_rpow2", "c[rfactor]i[nsize]i[nns]i[qual]i[etype]i[pscrn]i[cshift]s[fwidth]i" \
		"[fheight]i[ep0]f[ep1]f[threads]i[opt]i[fapprox]i[csresize]b[mpeg2]b[logicalCores]b[MaxPhysCore]b" \
		"[SetAffinity]b[threads_rs]i[logicalCores_rs]b[MaxPhysCore_rs]b[SetAffinity_rs]b[sleep]b" \
		"[prefetch]i[range]i[ThreadLevel]i[ThreadLevel_rs]i", Create_nnedi3_rpow2, 0);

	env->AddFunction("PointResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_PointResize, 0);
	env->AddFunction("BilinearResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_BilinearResize, 0);
	env->AddFunction("BicubicResizeMT", "c[target_width]i[target_height]i[b]f[c]f[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_BicubicResize, 0);
	env->AddFunction("LanczosResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_LanczosResize, 0);
	env->AddFunction("Lanczos4ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_Lanczos4Resize, 0);
	env->AddFunction("BlackmanResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_BlackmanResize, 0);
	env->AddFunction("Spline16ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_Spline16Resize, 0);
	env->AddFunction("Spline36ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_Spline36Resize, 0);
	env->AddFunction("Spline64ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_Spline64Resize, 0);
	env->AddFunction("GaussResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[p]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_GaussianResize, 0);
	env->AddFunction("SincResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_SincResize, 0);
	env->AddFunction("SinPowResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[p]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[ThreadLevel]i",FilteredResizeMT::Create_SinPowerResize, 0);

// Desample functions

	env->AddFunction("DeBilinearResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeBilinearResize, 0);
	env->AddFunction("DeBicubicResizeMT", "c[target_width]i[target_height]i[b]f[c]f[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeBicubicResize, 0);
	env->AddFunction("DeLanczosResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeLanczosResize, 0);
	env->AddFunction("DeLanczos4ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeLanczos4Resize, 0);
	env->AddFunction("DeBlackmanResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeBlackmanResize, 0);
	env->AddFunction("DeSpline16ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeSpline16Resize, 0);
	env->AddFunction("DeSpline36ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeSpline36Resize, 0);
	env->AddFunction("DeSpline64ResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeSpline64Resize, 0);
	env->AddFunction("DeGaussResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[p]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeGaussianResize, 0);
	env->AddFunction("DeSincResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[taps]i[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeSincResize, 0);
	env->AddFunction("DeSinPowResizeMT", "c[target_width]i[target_height]i[src_left]f[src_top]f[src_width]f[src_height]f[p]f[threads]i" \
		"[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[range]i[accuracy]i[order]i[ThreadLevel]i",FilteredResizeMT::Create_DeSinPowerResize, 0);

  env->AddFunction("aWarpSharp2", "c[thresh]i[blur]i[type]i[depth]i[chroma]i[depthC]i[cplace]s[blurV]i[depthV]i[depthVC]i" \
	  "[blurC]i[blurVC]i[threshC]i[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", Create_aWarpSharp, (void*)0);
  env->AddFunction("aWarpSharp", "c[depth]f[blurlevel]i[thresh]f[cm]i[bm]i[show]b" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", Create_aWarpSharp, (void*)1);
  env->AddFunction("aSobel", "c[thresh]i[chroma]i[threshC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", Create_aWarpSharp, (void*)2);
  env->AddFunction("aBlur", "c[blur]i[type]i[chroma]i[blurV]i[blurC]i[blurVC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", Create_aWarpSharp, (void*)3);
  env->AddFunction("aWarp", "cc[depth]i[chroma]i[depthC]i[cplace]s[depthV]i[depthVC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", Create_aWarpSharp, (void*)4);
  env->AddFunction("aWarp4", "cc[depth]i[chroma]i[depthC]i[cplace]s[depthV]i[depthVC]i" \
	  "[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i", Create_aWarpSharp, (void*)5);

    env->AddFunction("ConvertYUVtoLinearRGB",
		"c[Color]i[OutputMode]i[HDRMode]i[HLGLb]f[HLGLw]f[HLGColor]i[OOTF]b[EOTF]b[fullrange]b[mpeg2c]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertYUVtoLinearRGB, 0);
    env->AddFunction("ConvertLinearRGBtoYUV",
		"c[Color]i[OutputMode]i[HDRMode]i[HLGLb]f[HLGLw]f[HLGColor]i[OOTF]b[EOTF]b[fullrange]b[mpeg2c]b[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertLinearRGBtoYUV, 0);

    env->AddFunction("ConvertYUVtoXYZ",
		"c[Color]i[OutputMode]i[HDRMode]i[HLGLb]f[HLGLw]f[HLGColor]i[OOTF]b[EOTF]b[fullrange]b[mpeg2c]b" \
		"[Rx]f[Ry]f[Gx]f[Gy]f[Bx]f[By]f[Wx]f[Wy]f[Crosstalk]f" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertYUVtoXYZ, 0);
    env->AddFunction("ConvertXYZtoYUV",
		"c[Color]i[OutputMode]i[HDRMode]i[HLGLb]f[HLGLw]f[HLGColor]i[OOTF]b[EOTF]b[fullrange]b[mpeg2c]b[fastmode]b" \
		"[Rx]f[Ry]f[Gx]f[Gy]f[Bx]f[By]f[Wx]f[Wy]f[pColor]i[pRx]f[pRy]f[pGx]f[pGy]f[pBx]f[pBy]f[pWx]f[pWy]f[Crosstalk]f" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertXYZtoYUV, 0);

    env->AddFunction("ConvertRGBtoXYZ",
		"c[Color]i[OutputMode]i[HDRMode]i[HLGLb]f[HLGLw]f[HLGColor]i[OOTF]b[EOTF]b[fastmode]b" \
		"[Rx]f[Ry]f[Gx]f[Gy]f[Bx]f[By]f[Wx]f[Wy]f[Crosstalk]f" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertRGBtoXYZ, 0);
    env->AddFunction("ConvertXYZtoRGB",
		"c[Color]i[OutputMode]i[HDRMode]i[HLGLb]f[HLGLw]f[HLGColor]i[OOTF]b[EOTF]b[fastmode]b" \
		"[Rx]f[Ry]f[Gx]f[Gy]f[Bx]f[By]f[Wx]f[Wy]f[pColor]i[pRx]f[pRy]f[pGx]f[pGy]f[pBx]f[pBy]f[pWx]f[pWy]f[Crosstalk]f" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertXYZtoRGB, 0);

    env->AddFunction("ConvertXYZ_Scale_HDRtoSDR",
		"c[Coeff_X]f[Coeff_Y]f[Coeff_Z]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b" \
		"[prefetch]i[ThreadLevel]i",Create_ConvertXYZ_Scale_HDRtoSDR, 0);
    env->AddFunction("ConvertXYZ_Scale_SDRtoHDR",
		"c[Coeff_X]f[Coeff_Y]f[Coeff_Z]f[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b" \
		"[prefetch]i[ThreadLevel]i",Create_ConvertXYZ_Scale_SDRtoHDR, 0);

    env->AddFunction("ConvertXYZ_Hable_HDRtoSDR",
		"c[exposure_X]f[whitescale_X]f[a_X]f[b_X]f[c_X]f[d_X]f[e_X]f[f_X]f" \
		"[exposure_Y]f[whitescale_Y]f[a_Y]f[b_Y]f[c_Y]f[d_Y]f[e_Y]f[f_Y]f" \
		"[exposure_Z]f[whitescale_Z]f[a_Z]f[b_Z]f[c_Z]f[d_Z]f[e_Z]f[f_Z]f" \
		"[pColor]i[pRx]f[pRy]f[pGx]f[pGy]f[pBx]f[pBy]f[pWx]f[pWy]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertXYZ_Hable_HDRtoSDR, 0);
    env->AddFunction("ConvertRGB_Hable_HDRtoSDR",
		"c[exposure_R]f[whitescale_R]f[a_R]f[b_R]f[c_R]f[d_R]f[e_R]f[f_R]f" \
		"[exposure_G]f[whitescale_G]f[a_G]f[b_G]f[c_G]f[d_G]f[e_G]f[f_G]f" \
		"[exposure_B]f[whitescale_B]f[a_B]f[b_B]f[c_B]f[d_B]f[e_B]f[f_B]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertRGB_Hable_HDRtoSDR, 0);

    env->AddFunction("ConvertXYZ_Mobius_HDRtoSDR",
		"c[exposure_X]f[transition_X]f[peak_X]f[exposure_Y]f[transition_Y]f[peak_Y]f" \
		"[exposure_Z]f[transition_Z]f[peak_Z]f" \
		"[pColor]i[pRx]f[pRy]f[pGx]f[pGy]f[pBx]f[pBy]f[pWx]f[pWy]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertXYZ_Mobius_HDRtoSDR, 0);
    env->AddFunction("ConvertRGB_Mobius_HDRtoSDR",
		"c[exposure_R]f[transition_R]f[peak_R]f[exposure_G]f[transition_G]f[peak_G]f" \
		"[exposure_B]f[transition_B]f[peak_B]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertRGB_Mobius_HDRtoSDR, 0);

    env->AddFunction("ConvertXYZ_Reinhard_HDRtoSDR",
		"c[exposure_X]f[contrast_X]f[peak_X]f[exposure_Y]f[contrast_Y]f[peak_Y]f" \
		"[exposure_Z]f[contrast_Z]f[peak_Z]f" \
		"[pColor]i[pRx]f[pRy]f[pGx]f[pGy]f[pBx]f[pBy]f[pWx]f[pWy]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertXYZ_Reinhard_HDRtoSDR, 0);
    env->AddFunction("ConvertRGB_Reinhard_HDRtoSDR",
		"c[exposure_R]f[contrast_R]f[peak_R]f[exposure_G]f[contrast_G]f[peak_G]f" \
		"[exposure_B]f[contrast_B]f[peak_B]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertRGB_Reinhard_HDRtoSDR, 0);
		
    env->AddFunction("ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR",
		"c[Lhdr]f[Lsdr]f[CoeffAdj]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConvertLinearRGBtoYUV_BT2446_A_HDRtoSDR, 0);

    env->AddFunction("ConverXYZ_BT2446_C_HDRtoSDR",
		"c[ChromaC]b[PQMode]b[Lhdr]f[Lsdr]f[pct_ref]f[pct_ip]f[pct_wp]f[pct_sdr_skin]f[pct_hdr_skin]f[WhiteShift]f" \
		"[pColor]i[pRx]f[pRy]f[pGx]f[pGy]f[pBx]f[pBy]f[pWx]f[pWy]f[fastmode]b" \
		"[threads]i[logicalCores]b[MaxPhysCore]b[SetAffinity]b[sleep]b[prefetch]i[ThreadLevel]i",
		Create_ConverXYZ_BT2446_C_HDRtoSDR, 0);

	return PLUGINS_JPSDR_VERSION;	
}

