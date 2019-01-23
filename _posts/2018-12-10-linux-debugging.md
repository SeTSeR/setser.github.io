---
layout: post
author: Sergey
title: "Отладка программ в Linux"
comments: true
---

Довольно часто приходится сталкиваться с тем, что программа работает не так, как надо. В этой заметке я постараюсь описать наиболее простые инструменты и способы поиска ошибок в программе.

## GCC
Рассмотрим некоторые инструменты поиска ошибок, предоставляемые компилятором:
1. Предупреждения. GCC позволяет отлавливать большое количество ошибок посредством статического анализа, большая часть диагностик включается флагами компиляции -Wall -Wextra -Wpedantic.
2. Sanitizers. GCC умеет внедрять на этапе компиляции дополнительный инструментальный код, позволяющий отлавливать различные ошибки, например data races:
  * Address sanitizer. Позволяет отслеживать ошибки работы с памятью, такие как buffer-overflow, double free и use-after-free.
  * Thread sanitizer. Позволяет отслеживать различные ошибки, ведущие к race condition, например, вызов async-unsafe-функций в обработчике сигналов.
  * Leak sanitizer. Позволяет отслеживать утечки памяти.
  Весьма удобны в комбинации с отладочной информацией. При поиске ошибок по работе с памятью бывает полезен также valgrind.

## GDB
GDB - полноценный консольный отладчик для UNIX-like систем. Он позволяет останавливать выполнение в произвольных точках, отслеживать состояние переменных, стека, памяти и регистров и т. д. Возможна удалённая отладка.

Запуск:
```
gdb executable
```
Или просто:
```
gdb
```
После запуска можно указать файл для отладки:
```
file executable
```
Либо подключиться к удалённому сеансу отладки:
```
target remote host:port
```
Также можно отлаживать уже запущенную программу(нужны права суперпользователя):
```
attach <pid>
```
Запустить отлаживаемую программу на выполнение:
```
run arguments
```
Программа остановится, если выполнение дойдёт до точки останова, либо если программе придёт сигнал. Точки останова можно установить командой:
```
break main.c:15
```
либо
```
break 15
```
если отлаживается, например, программа из одного файла. Команды *step* и *next* выполняют следующую строку после остановки. Первая входит в функцию, если следующая строка содержит её вызов, вторая игнорирует вызов функции и останавливается уже после её выполнения. Команда *continue* продолжает выполнение до следующей точки останова, либо до следующего сигнала, либо до завершения программы. Команда *quit* осуществляет выход из GDB. С помощью команды
```
print выражение
```
можно выполнить любое выражение языка C, содержащее переменные, видимые в текущей точке программы. Команда *list* выдаёт листинг части программы в окрестности текущей строки.

GDB весьма удобен при отладке, если есть доступ к исходникам и отладочной информации. Если отладочной информации нет, для упрощения отладки можно использовать *objdump*, однако отладка всё равно будет затруднена вследствие оптимизаций компилятора. Возможно будет более удобно использование средств, предоставляемых ОС, таких как strace.

Что делать, если нет отладочной информации? GDB также поддерживает возможность отладки на уровне ассемблерных инструкций. Можно распечатать дизассемблированный листинг текущей функции с помощью команды *disassemble*. По умолчанию листинг выдаётся в синтаксисе AT & T, с помощью команды *set disassembly-flavor intel* это можно изменить. Навигация по инструкциям осуществляется с помощью команд *nexti*/*stepi*. Можно так же ставить точки останова по адресу и печатать значение регистров/памяти/информацию о текущем фрейме(*info frame*).

Для более подробной информации можно обратиться к документации GDB:
```
info gdb
```

## strace
strace - утилита, позволяющая отслеживать системные вызовы, используемые программой, и получаемые ей сигналы.
Запуск:
```
strace program
```
Либо:
```
strace -p pid
```
Для отладки уже запущенной программы(нужны права суперпользователя). Можно отлаживать многопроцессные программы, передав ключ *-f*. В этом случае будут печататься системные вызовы и сигналы не только родителя, но и детей, порождённых после присоединения к процессу. Можно указать файл, в который strace будет писать вывод:
```
strace -o log.txt program
```
Для системных вызовов strace печатает имя вызова, переданные ему аргументы и значение, которое он вернул. Можно ограничить количество отслеживаемых системных вызовов:
```
strace -e trace=open,close,read,write,fork,exec program
```
То же для сигналов:
```
strace -e signal=!SIGTERM program
```
Можно отслеживать данные, читаемые или записываемые в файловые дескрипторы:
```
strace -e read=0,3 program
strace -e write=1,2,4 program
```
Кроме того, можно вмешиваться в работу процесса, например, подменяя результаты системных вызовов. Об этом и многом другом можно прочитать в документации strace:
```
man 1 strace
```

## proc pseudo-filesystem
Различную отладочную информацию можно также почерпнуть из procfs. В Linux эта файловая система монтируется в /proc. Сведения о процессе с id pid можно найти в каталоге */proc/pid*.
Из такой информации можно выделить:
* Файл /proc/pid/exe является символической ссылкой на бинарный файл работающего процесса.
* Файл /proc/pid/status, содержащий различные данные о состоянии процесса, такие как:
    + Имя исполняемого файла.
    + umask процесса.
    + Собственно состояние процесса.
    + uid и gid процесса.
    + Количество потоков процесса.
    + Маски заблокированных, принимаемых и игнорируемых сигналов. Эти маски представлены в виде 8 байт, среди которых бит с номером i отвечает сигналу с номером i. Биты отсчитываются с единицы, справа налево.
* Файл /proc/pid/syscall содержит номера работающих в данный момент системных вызовов, вызванных из программы, и их аргументы.
* Файл /proc/pid/maps содержит список виртуальных страниц в адресном пространстве программы, файл /proc/pid/pagemap - их отображение в физические страницы.
* Файл /proc/pid/stack содержит трейс стека ядра.
* Каталог /proc/pid/fd содержит ссылки на файлы, соответствующие открытым в процессе файловым дескриптором, каталог /proc/pid/fdinfo содержит файлы с информацией об этих дескрипторах.
* Файл /proc/pid/limits содержит таблицу лимитов для процесса.

Более полное представление о структуре procfs можно получить в документации:
```
man 5 proc
```