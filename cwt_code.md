# CWT（连续小波变换）代码工程解析


- [CWT初始化](#cwt初始化)
  - [cwt\_init](#cwt_init)
  - [setCWTScales](#setcwtscales)
- [连续小波变换过程](#连续小波变换过程)
  - [步骤总览](#步骤总览)
  - [**Step1：初始化和参数设定**](#step1初始化和参数设定)
  - [**Step2：输入信号FFT**](#step2输入信号fft)


## CWT初始化

```c++
...

char *wave = "morlet"; // 使用Morlet小波，可选"paul"和"dog"
char *type = "pow";

cwt_object wt;

// 初始化参数
N = 504;
param = 6.0;
subscale = 2;
dt = 0.25;
s0 = dt;
dj = 1.0 / (double)subscale;
J = 11 * subscale;   // 总的尺度数
a0 = 2;            // 幂次

// 初始化小波变换对象
wt = cwt_init(wave, param, N, dt, J);

// 设置小波变换的尺度
setCWTScales(wt, s0, dj, type, a0);
```

其中`cwt_object`结构体定义如下：

<details>

```c++
    struct cwt_set{
        char wave[10];// Wavelet - morl/morlet,paul,dog/dgauss
        int siglength;// Length of Input Data
        int J;// Total Number of Scales
        double s0;// Smallest scale. It depends on the sampling rate. s0 <= 2 * dt for most wavelets
        double dt;// Sampling Rate
        double dj;// Separation between scales. eg., scale = s0 * 2 ^ ( [0:N-1] *dj ) or scale = s0 *[0:N-1] * dj
        char type[10];// Scale Type - Power or Linear
        int pow;// Base of Power in case type = pow. Typical value is pow = 2
        int sflag;
        int pflag;
        int npad;
        int mother;
        double m;// Wavelet parameter param
        double smean;// Input Signal mean

        cplx_data *output;
        double *scale;
        double *period;
        double *coi;
        double params[0];
    };
```

</details>


CWT结构体详细解析：

<details>

`cwt_set` 结构体旨在为连续小波变换（Continuous Wavelet Transform, CWT）配置参数并存储结果。以下是该结构体中各个字段的详细描述：

1. **wave[10]**：存储母小波函数的名称，如"morl"或"morlet"（莫雷特小波）、"paul"（保罗小波）、"dog"或"dgauss"（导数高斯小波）。

2. **siglength**：输入数据的长度。

3. **J**：小波变换使用的总尺度数。

4. **s0**：最小尺度，这个参数依赖于采样率。对于大多数小波，s0 应小于或等于 2 * dt。

5. **dt**：采样率，即数据采样的时间间隔。

6. **dj**：尺度之间的间隔。例如，尺度可以是 $ s0 \times 2^{[0:J-1] \times dj} $（对于指数类型的尺度变化）或 $ s0 \times [0:J-1] \times dj $（对于线性类型的尺度变化）。

7. **type[10]**：尺度类型，可以是"Power"或"Linear"，表示尺度变化是指数形式还是线性形式。

8. **pow**：如果尺度类型是"Power"，则此字段表示幂的基数，常见的值为2。

9. **sflag**：可能用于标记某种特定的状态或选项，具体用途可能需要根据上下文或其他文档来确定。

10. **pflag**：同样可能用于标记某种特定的状态或选项，具体用途可能需要查阅相关文档。

11. **npad**：填充后的数据长度，通常用于FFT计算中，以提高频率分辨率。

12. **mother**：用于标识选用的小波母函数的编号或类型。

13. **m**：小波参数，用于调整小波函数的形状，具体含义取决于所选的小波母函数。

14. **smean**：输入信号的均值，用于在处理前去除信号的直流分量。

15. **output**：指向复数数据类型（cplx_data）的指针，用于存储小波变换的输出结果。

16. **scale**：存储每个尺度的实际值的数组。

17. **period**：存储每个尺度对应的周期的数组。

18. **coi**：存储每个时间点的锥形影响区（Cone of Influence）的数组。

19. **params[0]**：一个大小可变的数组，用于存储可能需要的额外参数，其大小和内容依赖于具体应用。

该结构体为执行小波变换提供了必要的配置信息，并存储了变换的结果，以便于进一步的分析和处理。

</details>

### cwt_init
[原始代码](./src/wavelib.c#L283)

初始化一个小波变换对象 cwt_object 的函数 cwt_init。这个函数根据指定的小波类型、参数、信号长度、采样率和尺度数来配置小波变换的各种参数。

<details>

```c
/**
 * @description: 连续小波变换对象初始化
 * @param {char*} wave: 指定小波函数的名称（如 "morlet", "paul", "dog"）。
 * @param {double} param: 特定小波函数的参数。
 * @param {int} siglength: 输入信号的长度。
 * @param {double} dt: 采样间隔。
 * @param {int} J: 小波变换使用的总尺度数。
 * @return {*}
 */
cwt_object cwt_init(const char* wave, double param,int siglength, double dt, int J) {
	cwt_object obj = NULL;
	int N, i,nj2,ibase2,mother;
	double s0, dj;
	double t1;
	int m, odd;
	const char *pdefault = "pow";

	m = (int)param;
	odd = 1;
	if (2 * (m / 2) == m) {
		odd = 0;
	}

	N = siglength;
	nj2 = 2 * N * J;

	// 构造cwt对象
	obj = (cwt_object)malloc(sizeof(struct cwt_set) + sizeof(double)* (nj2 + 2 * J + N));

	if (!strcmp(wave, "morlet") || !strcmp(wave, "morl")) {
		s0 = 2 * dt;
		dj = 0.4875; // why?
		mother = 0;
		if (param < 0.0) {
			printf("\n Morlet Wavelet Parameter should be >= 0 \n");
			exit(-1);
		}
		if (param == 0) {
			param = 6.0;
		}
		strcpy(obj->wave,"morlet");
		
	}
	else if (!strcmp(wave, "paul")) {
		s0 = 2 * dt;
		dj = 0.4875;
		mother = 1;
		if (param < 0 || param > 20) {
			printf("\n Paul Wavelet Parameter should be > 0 and <= 20 \n");
			exit(-1);
		}
		if (param == 0) {
			param = 4.0;
		}
		strcpy(obj->wave,"paul");
	
	}
	else if (!strcmp(wave, "dgauss") || !strcmp(wave, "dog")) {
		s0 = 2 * dt;
		dj = 0.4875;
		mother = 2;
		if (param < 0 || odd == 1) {
			printf("\n DOG Wavelet Parameter should be > 0 and even \n");
			exit(-1);
		}
		if (param == 0) {
			param = 2.0;
		}
		strcpy(obj->wave,"dog");
	}

	obj->pow = 2;
	strcpy(obj->type, pdefault);

	obj->s0 = s0; 	// 	最小尺度
	obj->dj = dj;	// 	尺度之间的间隔
	obj->dt = dt;	// 	采样率，即数据采样的时间间隔
	obj->J = J;		// 	小波变换使用的总尺度数
	obj->siglength = siglength; 	// 	输入信号长度
	obj->sflag = 0;			//
	obj->pflag = 1;			//
	obj->mother = mother;	// 	母小波
	obj->m = param;			//	小波参数，用于调整小波函数的形状

	t1 = 0.499999 + log((double)N) / log(2.0);
	ibase2 = 1 + (int)t1;

	obj->npad = (int)pow(2.0, (double)ibase2);

	// 分配数据地址
	obj->output = (cplx_data*) &obj->params[0];
	obj->scale = &obj->params[nj2];
	obj->period = &obj->params[nj2 + J];
	obj->coi = &obj->params[nj2 + 2*J];

	for (i = 0; i < nj2 + 2 * J + N; ++i) {
		obj->params[i] = 0.0;
	}

	return obj;
}

```

</details>

### setCWTScales

设置小波尺度等参数。

<details>

```c
void setCWTScales(cwt_object wt, double s0, double dj,const char *type,int power) {
	int i;
	strcpy(wt->type,type);
	//s0*pow(2.0, (double)(j - 1)*dj);
	if (!strcmp(wt->type, "pow") || !strcmp(wt->type, "power")) {
		for (i = 0; i < wt->J; ++i) {
			wt->scale[i] = s0*pow((double) power, (double)(i)*dj);
		}
		wt->sflag = 1;
		wt->pow = power;
		
	}
	else if (!strcmp(wt->type, "lin") || !strcmp(wt->type, "linear")) {
		for (i = 0; i < wt->J; ++i) {
			wt->scale[i] = s0 + (double)i * dj;
		}
		wt->sflag = 1;
	}
	else {
		printf("\n Type accepts only two values : pow and lin\n");
		exit(-1);
	}
	wt->s0 = s0;
	wt->dj = dj;
}
```

</details>


## 连续小波变换过程

[原始代码](./src/wavelib.c#L1576)

跳过执行小波变换前的一些预处理步骤，详细解析[小波变换主体](./src/cwt.c#L170)过程。

### 步骤总览

- **初始化和参数设定**

   1. **函数签名**：cwavelet函数接受原始数据、时间序列长度、时间间隔、小波母函数类型、小波参数、最小尺度、尺度间隔、总尺度数、填充长度、输出的小波变换结果、尺度数组、周期数组和锥形影响区（Cone of Influence, COI）数组。

   2. **参数检查**：检查npad（填充后的长度）是否大于或等于N（原始数据长度），若不满足则打印错误信息并退出。
   
   3. **FFT对象初始化**：初始化两个FFT对象，一个用于正向FFT，另一个用于逆向FFT。
   
   4. **内存分配**：为填充后的数据、FFT结果、小波“女儿”波形和波数数组分配内存。

- **数据预处理**

    5. **计算输入数据的均值**：计算原始数据的均值，用于去除直流分量。

    6. **填充和去直流**：用原始数据减去均值填充ypad数组的前N个元素，剩余部分填充0，以准备进行FFT。

- **FFT变换**

    7. **执行FFT**：对去直流后的数据执行FFT，得到频域表示。

    8. **归一化FFT结果**：将FFT的结果除以npad以进行归一化。

    9. **构建波数数组**：计算每个频率点对应的波数，用于后续的小波变换。

- **小波变换主循环**

    10. **遍历所有尺度**：对每一个尺度，使用小波母函数和对应的波数计算小波变换。

    11. **小波函数生成**：根据当前尺度和波数，计算小波“女儿”波形。

    12. **小波变换**：将小波“女儿”波形与FFT的结果相乘，得到小波变换的频域表示。

    13. **执行逆FFT**：对小波变换的结果执行逆FFT，得到时域的小波变换结果。

- **结果处理和内存释放**

    14. **记录结果**：将每个尺度的小波变换结果保存到输出数组中。

    15. **计算锥形影响区**：根据小波变换的特性，计算并记录每个时间点的锥形影响区。

    16. **释放内存**：释放所有分配的内存资源。

    17. **清理FFT对象**：释放FFT对象。




### **Step1：初始化和参数设定**

```c
  int i, j, k, iter;
  double ymean, freq1, pi, period1, coi1;
  double tmp1, tmp2;
  double scale1;
  double *kwave;
  fft_object obj, iobj;
  fft_data *ypad, *yfft, *daughter;

  (void)s0;
  (void)dj; /* yes, we need these parameters unused */

  pi = 4.0 * atan(1.0);

  if (npad < N) {
    printf("npad must be >= N \n");
    exit(-1);
  }

  obj = fft_init(npad, 1);
  iobj = fft_init(npad, -1);


  ypad = (fft_data *)malloc(sizeof(fft_data) * npad);
  yfft = (fft_data *)malloc(sizeof(fft_data) * npad);
  daughter = (fft_data *)malloc(sizeof(fft_data) * npad);
  kwave = (double *)malloc(sizeof(double) * npad);
```

### **Step2：输入信号FFT** 
  
```c
  // 以下是FFT标准流程：计算均值->去直流->执行FFT
  ymean = 0.0;

  for (i = 0; i < N; ++i) {
    ymean += y[i];
  }

  ymean /= N;

  for (i = 0; i < N; ++i) {
    ypad[i].re = y[i] - ymean;
    ypad[i].im = 0.0;
  }

  for (i = N; i < npad; ++i) {
    ypad[i].re = ypad[i].im = 0.0;
  }

  // Find FFT of the input y (ypad)
  fft_exec(obj, ypad, yfft);
  for (i = 0; i < npad; ++i) {
    yfft[i].re /= (double)npad;
    yfft[i].im /= (double)npad;
  }
```

