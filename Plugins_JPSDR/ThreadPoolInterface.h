#ifndef __ThreadPoolInterface_H__
#define __ThreadPoolInterface_H__

#include <Windows.h>

#include "ThreadPoolDef.h"

#define THREADPOOLINTERFACE_VERSION "ThreadPoolInterface 1.6.1"

typedef struct _UserData
{
	volatile uint16_t UserId;
	volatile int8_t nPool;
} UserData;


class ThreadPoolInterface
{
	public :

	virtual ~ThreadPoolInterface(void);
	static ThreadPoolInterface* Init(uint8_t num);

	uint8_t GetThreadNumber(uint8_t thread_number,bool logical);
	int16_t AddPool(uint8_t num);
	bool CreatePool(uint8_t num);
	bool DeletePool(uint8_t num);
	bool RemovePool(uint8_t num);
	bool AllocateThreads(uint16_t &UserId,uint8_t thread_number,uint8_t offset_core,uint8_t offset_ht,bool UseMaxPhysCore,bool SetAffinity,bool sleep,int8_t nPool);
	bool ChangeThreadsAffinity(uint16_t UserId,uint8_t offset_core,uint8_t offset_ht,bool UseMaxPhysCore,bool SetAffinity,bool sleep,int8_t nPool);
	bool DeAllocateThreads(uint16_t UserId);
	bool RequestThreadPool(uint16_t UserId,uint8_t thread_number,Public_MT_Data_Thread *Data,int8_t nPool,bool Exclusive);
	bool ReleaseThreadPool(uint16_t UserId,bool sleep);
	bool StartThreads(uint16_t UserId);
	bool WaitThreadsEnd(uint16_t UserId);
	bool GetThreadPoolStatus(uint16_t UserId,int8_t nPool);
	uint8_t GetCurrentThreadAllocated(uint16_t UserId,int8_t nPool);
	uint8_t GetCurrentThreadUsed(uint16_t UserId,int8_t nPool);
	uint8_t GetLogicalCPUNumber(void);
	uint8_t GetPhysicalCoreNumber(void);

	protected :
	
	volatile bool Status_Ok;
	volatile uint8_t NbrePool;

	public :

	bool GetThreadPoolInterfaceStatus(void) {return(Status_Ok);}
	int8_t GetCurrentPoolCreated(void) {return((Status_Ok) ? NbrePool:-1);}

	protected :

	ThreadPoolInterface(void);

	CRITICAL_SECTION CriticalSection;
	HANDLE ghMutexResources;
	BOOL CSectionOk;
	HANDLE JobsEnded[MAX_THREAD_POOL],ThreadPoolFree[MAX_THREAD_POOL];
	UserData TabId[MAX_USERS];
	volatile uint16_t NbreUsers;
	HANDLE EndExclusive;
	volatile bool Error_Occured;

	volatile bool ThreadPoolRequested[MAX_THREAD_POOL],JobsRunning[MAX_THREAD_POOL];
	volatile bool ThreadPoolReleased[MAX_THREAD_POOL],ThreadWaitEnd[MAX_THREAD_POOL];
	volatile bool ThreadPoolWaitFree[MAX_THREAD_POOL];
	volatile bool ExclusiveMode;
	volatile uint8_t NbrePoolEvent;

	bool CreatePoolEvent(uint8_t num);
	void FreeData(void);
	void FreePool(void);
	bool EnterCS(void);
	void LeaveCS(void);
	bool GetMutex(void);
	void FreeMutex(void);
	
	private :

	ThreadPoolInterface (const ThreadPoolInterface &other);
	ThreadPoolInterface& operator = (const ThreadPoolInterface &other);
	bool operator == (const ThreadPoolInterface &other) const;
	bool operator != (const ThreadPoolInterface &other) const;
};

#endif // __ThreadPoolInterface_H__

