# Matlab scripts for TOPPE v2b

These scripts work with the toppev2b.e driver/psd.

These scripts are self-contained; the only dependency is the ./+utils/ folder.

For conversion to/from Pulseq, use ge2seq.m/seq2ge.m

https://toppemri.github.io/


## A typical workflow

1. Design RF and gradient waveforms (functions in +utils/+rf/ and +utils/+epi/ may be helpful for you)

2. Edit systemspecs.m according to your scanner, or call with options.
```
>> sys = toppe.systemspecs();
```

3. Use 'writemod.m' to write each set of rf/gradient waveforms to a .mod file, e.g.,
```
>> toppe.writemod('rf', rfwav, 'gz', gzwav, 'ofname', 'tipdown.mod', 'system', sys);
```

4. Create 'scanloop.txt' with detailed scan loop instructions.
   See 'writeloop.m' in any of the examples folders.

5. Copy modules.txt, scanloop.txt, and the .mod files to /usr/g/bin/ on the scanner, and scan with toppev2.e

Many complete examples are provided in matlab/examples


## More usage examples:

Write a slice-selective RF excitation module on a system with max slew rate 130 T/m/s:
```
>> sys = toppe.systemspecs('MaxSlew', 130, 'slewUnit', 'T/m/s');
>> toppe.writemod('rf', rfwav, 'gz', gzwav, 'ofname', 'tipdown.mod', 'system', sys);
```

Display part of the sequence:
```
>> ibeg = 10;
>> iend = 20;
>> toppe.plotseq(ibeg, iend);      % by default, looks for 'modules.txt' and 'scanloop.txt' in current folder
>> d = toppe.readloop('scanloop.txt');
>> sys = toppe.systemspecs();
>> sys.toppe.timessi = 200;   % us
>> toppe.plotseq(ibeg, iend, 'loopArr', d, 'system', sys);
```

Display sequence in movie/loop mode.
```
>> nModulesPerTR = 3;
>> toppe.playseq(nModulesPerTR);
```

## Using +toppe package functions in your .m files

Add the *root* folder (the folder containing the +toppe folder) to your Matlab path and either:

1. Call functions with "toppe." prefix, e.g.,
```
[rf,gx] = toppe.readmod('module1.mod');
```

or

2. Add the following to the top of your .m file:
```
import toppe.*
import toppe.utils.*
```
and call the functions without the prefix:
```
[rf,gx] = readmod('module1.mod');
```

To list package contents, do
```
>> help toppe
```

## Changes from toppev1

* It's now possible to associate multiple waveforms with each .mod file; accordingly, an additional column is required in scanloop.txt specifying which of the waveform(s) to play out.

* The header in the .mod file has changed, so be sure to use the writemod.m in this folder when scanning with toppev2.

* Changed the interface/usage of most of the .m functions; now uses the popular ('param', value) syntax.

* writemod.m now includes checks for system hardware limits, and should be more robust.

* getscantime.m: returns total scan time (new)

## Variable and function naming convention

* function names are all lower case

* variable and struct names are mixed case starting with lower case (e.g., maxLength)

* The name of an array of variables/structs is <variable name>Array (e.g., moduleArr)

* prefix n is used for variables representing the number of objects

