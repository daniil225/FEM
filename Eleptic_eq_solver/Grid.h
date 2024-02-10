#pragma once

#include <string>

/* GridBuilder v1.0.0 */

/* Возникла идея всю эту систему в следующей версии определить в класс 
	Причин для этого предостаточно 
	1) Явно идет разделение на приватный и публичный функционал 
	2) Возможность организовать паттерн управления памятью RAII
	3) Внутреннее представление структуры сетки никак не должно влиять на конечного пользовотеля 
	4) Большое количество специфических струтур неоюходимых для ораганизации структуры взаимоедйствия можно определить всего несколько для доступа к внешним данным 
	5) Плюс явно проглядывается базовая структура для интерфеса 
*/


/* Чтение и сетки из файла. Разбиение сетки.  */
/* Представленная структура используется для хранения сетки в виде шестигранных элементов 
   Прямые призмы с четырехугольным основанием с основанием в плоскости (x;z)
*/

/* Формат файла 
* Nx Nz
* x1,1 z1,1 x2,1 z1,2 ... xn,1 z1,n
* x1,2 z2,1 x2,2 z2,2 ... x2,n z2,n
* .....
* Ny
* y1, y2, ..., yk
* 
* L - целое число - количество подобластей и далее L наборов чисел по 7 штук 
* 1 число - номер формул определяющих параметры ДУ в подобласти
* 2 число - первая вертикальная ломанная определеяет правую границу 
* 3 число - вторая вертикальная ломанная определяет левую границу 
  их номера определяются в соответствии с правилом для k узла их номера в массивах точек будут: k, (Nx + k), (2*Nx + k), ... , ( (Ny-1)*Nx + k )
* 4 число - опеределяет горизонтальную границу снизу
* 5 число - определяет горизонтальную границу сверху 
* 6 число - граница по y начало
* 7 число - граница по y конец
* Дальше 3 сторики задающие разбиение сетки 
* P - целое число - количество подоблостей граничных условий и далее P наборов чисел по 8 штук
* 1 число - номер формул определяющих параметры граничных условий 
* 2 число - Тип краевого условия 0 - соответсвует не заданным краевым 
* 3 число - первая вертикальная ломанная определеяет правую границу 
* 4 число - вторая вертикальная ломанная определяет левую границу 
* 5 число - опеределяет горизонтальную границу снизу
* 6 число - определяет горизонтальную границу сверху 
* 7 число - граница по y начало
* 8 число - граница по y конец
* !!! Отличительная особенность одна из координатных линий должна быть фиксированной 
* N1 coef1 N2 coef 2 ...
* ...
* ...
* Первый набор задает разбиение каждого отрезка по х
* Второй набора задает разбиение по z
* Третий набор задает разбиение по у
*/

namespace Grid
{

	/* Секция с функциями и структурами определяющими дополнительную информацию о Точке в области  */
	namespace PointInfo
	{
		/*
		Структура определяющая информацию о точке в области
		Сжатая структура
		*/
		struct Info
		{
			uint64_t BaseInfo : 8; /* Базовая информация */
			uint64_t AreaInfo : 32; /* Информация об области определяющей параметры ДУ */
			uint64_t BoundInfo : 16; /*Краевые условия просто набор формул */
			uint64_t TypeBoundCond : 8; /* Тип Краевого условия */
		};

		/* Развертка информации о типе КУ  */
		struct BoundInfo
		{
			uint8_t size = 0; // Сколько типов задано 
			uint8_t Cond[4]{ 0,0,0,0 }; // Набор формул определяющих КУ
			uint8_t TypeCond[4]{ 0,0,0,0 }; // Тип Краевого условия 
		};

		/* Развертка информации об области */
		struct AreaInfo
		{
			uint8_t size = 0;
			uint8_t Cond[8]{ 0,0,0,0,0,0,0,0 };
		};

		struct Comand
		{
			typedef const uint32_t ComandType;
			/* Clear comand */
			static ComandType BaseInfoClear = 0x00;
			static ComandType AreaInfoClear = 0x00000000;
			static ComandType BoundInfoClear = 0x0000;
			static ComandType TypeBoundCondClear = 0x00;

			/* Base Info comand */
			static ComandType SetZeroFiFictitious = 0xFE;
			static ComandType GetFiFictitious = 0x01;

			/* Area Comand */
			static ComandType SetZeroAreaCount = 0xE1;
			static ComandType GetAreaCount = 0x1E;
			static ComandType ShiftAreaCountBits = 1;
			static ComandType RightBoundAreaCount = 8;
			static ComandType RightBoundValAreaInfo = 15;
			static ComandType BaseAreaShift = 4;
			static ComandType GetAreaInfoBaseComand = 0x0000000F;

			/* Bound Comand */
			static ComandType SetZeroBountCount = 0x1F;
			static ComandType ShiftBoundCountBits = 5;
			static ComandType GetBoundCount = 0xE0;
			static ComandType RightBoundBoundCount = 4;
			static ComandType RightBoundValBoundInfo = 15;
			static ComandType RightBoundValBoundType = 3;
			static ComandType BaseBoundShift = 4;
			static ComandType BaseTypeBoundShift = 2;
			static ComandType GetBoundInfoBaseComand = 0x000F;
			static ComandType GetTypeBoundBaseComand = 0x03;
		};

		/* Функции управления битовыми структурами */
		/* Private Section */

		/* Очистка структуры */
		void ClearInfo(Info& info);

		/* Установка фиктивный/не фиктивный узел 
			@param
			Info& info - информационная структура 
			uint8_t val true - не фиктивный false - фиктивный
			ret void
		*/
		void SetFictitiousBit(Info& info, uint8_t val);
		
		/* Установка Количества различныйх типов областей которым принадлежит точка 
			@param
			Info& info - информационная структура 
			uint8_t val - диапазон [0;8]
			ret void
		*/
		void SetAreaCountBits(Info& info, uint32_t val);

		/* Установка Количества различных типов границ к которым примыкает данная точка
			@param
			Info& info - информационная структура
			uint8_t val - диапазон [0;4]
			ret void
		*/
		void SetBoundCountBits(Info& info, uint32_t val);

		/* Получить количество различных областей к которым примыкает точка
			@param
			Info& info - информационная структура
			ret uint32_t - информациоя о том к скольки областям примыкает  диапазон [0;8]
			0 - ни к одной 
			1 - к одной 
			...
		*/
		uint32_t GetAreaCount(const Info& info);

		/* Получить количество различных областей к которым примыкает граница 
			@param
			Info& info - информационная структура
			ret uint32_t - информациоя о том к скольки Границам примыкает(различным)  диапазон [0;4]
			0 - ни к какой 
			1 - к одной 
			2 - 
			3 - 
			4 - 
		*/
		uint32_t GetBoundCount(const Info& info);
		/*********************************************/

		/* Public section */

		/*
			@param
			Info& info 
			uint32_t val - набор формул задающий значения на границе диапазон [0,15]
			uint32_t boundType - тип Ку диапазон  [0,3]
		*/
		void SetBoundInfo(Info& info, uint32_t val, uint32_t boundType);
		/*
			@param
			Info& info
			uint32_t val - набор формул задающий значения в области определяющей параметры ДУ  диапазон [0,15]
		*/
		void SetAreaInfo(Info& info, uint32_t val);

		/* 
			@param
			const Info& info
			ret bool - false - фиктивный true - не фиктивный 
		*/
		bool IsFiFictitious(const Info& info);

		/* 
			@param
			const Info& info
			ret bool - false - не граница true - граница 
		*/
		bool IsBound(const Info& info);

		/*
			@param
			const Info& info
			ret BoundInfo -структура содержащая информацию о типах Ку и формул которые их задают
		*/
		BoundInfo GetBoundInfo(const Info& info);

		/*
			@param
			const Info& info
			ret AreaInfo -структура содержащая информацию о типах областей которым принадлежит точка и как следствие это определяет парметры ДУ 
		*/
		AreaInfo GetAreaInfo(const Info& info);
		/********************************************/

		/* Debug functions */
		void PrintInfo(const Info& info);

		void PrintBoundInfo(const BoundInfo& Bound);

		void PrintAreaInfo(const AreaInfo& Area);
	};
	/***********************************************************************************/

	/*  Структура данных для представления области  */
	const double eps = 1e-13; // Машинный ноль 

	/* Струкутура  базовой сетки на ее основе генерируется основная сетка для расчета */
	struct BaseGrid
	{
		int CountOfDivision = 1; // Количество делений базовой настройки сетки; 1 - соответсвует тому что делений не было 2 тому что исходная область разделена в двое и.т.д
		const int SizeOfCalculationAreaElemet = 7;
		const int SizeOfBoundsAreaElement = 8;
		const int SizeOfDivideParam = 3;

		/* Структура для базовой точки  */
		struct PointXZS
		{
			double x;
			double z;
		};

		struct DivideParamS
		{
			int num; // количество интервалов на которое нужно разделить отрезок 
			double coef; // Коэффициент растяжения или сжатия 
		};
		int Nx = 0; // количество узлов вдоль горизонатального направления 
		int Ny = 0; // количество узлов вдоль вертикального направления 
		int Nz = 0; // количество узлов вдоль оси z
		int L = 0; // Количество подоблостей 
		int P = 0; // Количество видов границ 
		PointXZS **BaseGridXZ = nullptr; //  Базовая сетка в плосткости XZ 
		double *BaseGridY = nullptr; //  Базовая сетка по Y
		/*
			| номер формул | граница по x(start) | граница по x(end) | | граница по z(start) | граница по z(end) | граница по y(start) | граница по y(end) |
		*/
		int** CalculationArea = nullptr; // Массив расчетных областей
		/*
			| номер формул | | тип КУ | граница по x(start) | граница по x(end) | | граница по z(start) | граница по z(end) | граница по y(start) | граница по y(end) |
		*/
		int** BoundsArea = nullptr; // Границы

		/*
			3 массива 
			DivideParam[0] - массив для разбиения по оси x размер Nx-1
			DivideParam[1] - массив для разбиения по оси z размер Nz-1
			DivideParam[2] - массив для разбиения по оси y размер Ny-1
		*/
		DivideParamS **DivideParam = nullptr; // Массив для разбиения 
	};



	// Геометрическая точка или узел конечного элемента 
	struct Point
	{
		// Информация о границе конкретный набор формул + информация об узле фиктивный или нет + информация о том граничный ли узел или нет 
		// + информация о конечном элементе + тип КУ + точка в пространтве 
		// Вся эта информация хранится в соответсвующих битах числа битовая картат ниже
		
		PointInfo::Info info;

		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
	};


	// Конечный элемент 
	struct Element
	{
		struct BoundS
		{
			bool IsBound = false;
			int LocalIdx[4]{ -1,-1,-1,-1 };
			int GlobalIdx[4]{ -1, -1, -1, -1 };
			int BoundType = -1;
			int BoundInfo = -1;
		};
		// Локальные узлы на элементе в них лежат глобальные координаты, но для построения алгоритма обработки лучше возвращать такой елемент 
		// Из такой структуры можно легко получить всю информацию

		static const int FinitElementSize = 8; // Размер элемента 
		static const int BoundCount = 6;
		/* Информация о конечном элементе */
		bool isFictive = false; // по умолчанию фиктивный 
		int AreaInfo = -1; // Номер формулы задающий парметры ДУ на заданной области 
		/* Граница */
		/* Всего на конечном элементе может быть 6 границ (8-и угольник) 
			Со следующей локальной нумерацией 
			{1,2,3,4}, {5,6,7,8}, {1,3,5,7}, {2,4,6,8}, {1,2,5,6}, {3,4,7,8} - в нумерации с 1
			на каждой границе могут быть заданы свои условия. При этом краевое условие определяется 2 числоми 
			1 число - номер формул задающих КУ
			2 число - Тип КУ
			Заведем 6 структур для каждой границы отдельно 

			Структура в следующем формате:
			bool IsBound = false; - Граница ли это. По умолчанию нет 
			int LocalIdx[4]; - Локальная нумерация 
			int GlobalIdx[4]; - Соответсвующая глобальная нумерация для внесения результатов в Глобальную матрицу 
			int BoundType = -1; - Тип КУ
			int BoundInfo = -1; - Набор фомул 
		*/
		BoundS Bound[BoundCount]; // Вся информация о конечном элементе 

		Point e[FinitElementSize];
		int GlobalIdx[FinitElementSize];
	};

	struct GridS
	{
		int Dim = 0; // Размерность сетки и СЛАУ в то же время  
		int Nx = 0; // Сумарное количество узлов по оси Х
		int Nz = 0; // Сумарное количество узлов по оси Z
		int Ny = 0; // Сумарное количество узлов по оси У все величины получаются после генерации сетки 
		Point* Grid = nullptr; // Массив точек получающийся при генерации конечных элементов 
	};

	/* Загрузка данных из файла 
		@param
		std::string filename - имя файла хранящего структуру сетки  
		BaseGrid &baseGrid - Базовая сетка по которой будет генерироваться Основаная расчетная область 
		@ret: void
	*/
	void Load(std::string filename, BaseGrid &baseGrid);

	/*  Генерация расчетной области
		@param
		const BaseGrid &baseGrid - базовая сетка описания расчетной области 7
		Grid &grid - Сгенерированная расчетная область 
		ret: void
	*/
	void GenerateGrid(BaseGrid &baseGrid, GridS &grid);


	/* Перегенерация области, если она уже была задана
		@param
		const BaseGrid &baseGrid - базовая сетка описания расчетной области 7
		Grid &grid - Сгенерированная расчетная область
		ret: void
		!!! Происходит очистка памяти, если она будет не выделена будет ошибка 
	*/
	void ReGenerateGrid(BaseGrid& baseGrid, GridS& grid);


	/* Возвращает конечный элемент по его номеру 
		@param 
		const  GridS& grid - сетка 
		int idx
		ret Element - структура из 8-и точек с локальной нумерацией. Каждая точка содержит соответсвующую информацию о том какой границе принадлежит 
	*/
	Element GetFinitElement(const GridS& grid,int idx);

	/* Дробление сетки на заданный коэффициент 
	* !!! Вызывать строго до генерации области 
		@param
		BaseGrid& baseGrid - базовая сетка (в этом поле будем умнажать количество подинтервалов 
		int coef
		ret void
		результат раздробленная сетка в заданное количесво раз 
	*/
	void DivideGrid(BaseGrid& baseGrid, int coef);


	/* Очистка памяти структуры основной сетки
		@param
		GridS& grid - стека
		ret void
	*/
	void DeallocateGrid(GridS& grid);


	/* Деалокаторы памяти */
	void DeallocateBaseGrid(BaseGrid& baseGrid);

	/* Debug function */
	void PrintListGrid(const GridS& grid, int list);

	void PrintElement(const Element& element);

	/********************************************************************************************************/

	/* Приватная секция в ней происходит манипуляции с областью на низком уровне */

	namespace Grid_private
	{
		/* Расчитывает общее число узлов которое получится в КЭ сетке */
		/*
			@param
			const BaseGrid &baseGrid - 
			Grid &grid - 
			ret void 
			результат в grid поля Nx,Ny,Nz инициализируются своими значениями 
		*/
		void GetTotalNumberOfNodes(const BaseGrid &baseGrid, GridS &grid);
		
		/* Генерация всей расчетной области без учетка фиктивных элементов и принадлежности к какой либо границе и расчетной области */
		/*
			@param
			const BaseGrid &baseGrid -
			Grid &grid -
			ret void
			результат в grid поле Grid получит свои значения узлов 
		*/
		void GenerateBaseGrid(const BaseGrid& baseGrid, GridS& grid);

		/* Функция учетка принадлежности к границе и учетка Факта Фиктивный или нет узел */
		/*
			@param
			const BaseGrid &baseGrid -
			Grid &grid -
			ret void
			результат в grid поле Points.Info  получит информацию о том принадлежит ли точка какой-либо области  + номер формул для данной области 
		*/
		void DivisionIntoSubAreas(const BaseGrid& baseGrid, GridS& grid);

		/* Функция учетка типа КУ и установка факта является ли элемент граничным */
		/*
			@param
			const BaseGrid &baseGrid -
			Grid &grid -
			ret void
			результат в grid поле Points.Info  получит информацию о том принадлежит ли точка какой-либо границе + номер формул для данной области
		*/
		void DivisionIntoSubBounds(const BaseGrid& baseGrid, GridS& grid);


	};

};