у·
°▄
B
AssignVariableOp
resource
value"dtype"
dtypetypeИ
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
Ы
Conv2D

input"T
filter"T
output"T"
Ttype:	
2"
strides	list(int)"
use_cudnn_on_gpubool(",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 "-
data_formatstringNHWC:
NHWCNCHW" 
	dilations	list(int)

W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(И

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetypeИ
E
Relu
features"T
activations"T"
Ttype:
2	
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0И
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0И
?
Select
	condition

t"T
e"T
output"T"	
Ttype
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
H
ShardedFilename
basename	
shard

num_shards
filename
9
Softmax
logits"T
softmax"T"
Ttype:
2
N
Squeeze

input"T
output"T"	
Ttype"
squeeze_dims	list(int)
 (
╛
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring И
@
StaticRegexFullMatch	
input

output
"
patternstring
Ў
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
-
Tanh
x"T
y"T"
Ttype:

2
Ц
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 И"serve*2.4.12unknown8╓║
|
conv1d/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:░	а*
shared_nameconv1d/kernel
u
!conv1d/kernel/Read/ReadVariableOpReadVariableOpconv1d/kernel*$
_output_shapes
:░	а*
dtype0
o
conv1d/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:а*
shared_nameconv1d/bias
h
conv1d/bias/Read/ReadVariableOpReadVariableOpconv1d/bias*
_output_shapes	
:а*
dtype0
w
dense/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:аГР*
shared_namedense/kernel
p
 dense/kernel/Read/ReadVariableOpReadVariableOpdense/kernel*!
_output_shapes
:аГР*
dtype0
m

dense/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:Р*
shared_name
dense/bias
f
dense/bias/Read/ReadVariableOpReadVariableOp
dense/bias*
_output_shapes	
:Р*
dtype0
y
dense_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	Р*
shared_namedense_1/kernel
r
"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel*
_output_shapes
:	Р*
dtype0
p
dense_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_1/bias
i
 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
_output_shapes
:*
dtype0
x
dense_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:<*
shared_namedense_2/kernel
q
"dense_2/kernel/Read/ReadVariableOpReadVariableOpdense_2/kernel*
_output_shapes

:<*
dtype0
p
dense_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:<*
shared_namedense_2/bias
i
 dense_2/bias/Read/ReadVariableOpReadVariableOpdense_2/bias*
_output_shapes
:<*
dtype0
y
dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	<░	*
shared_namedense_3/kernel
r
"dense_3/kernel/Read/ReadVariableOpReadVariableOpdense_3/kernel*
_output_shapes
:	<░	*
dtype0
q
dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:░	*
shared_namedense_3/bias
j
 dense_3/bias/Read/ReadVariableOpReadVariableOpdense_3/bias*
_output_shapes	
:░	*
dtype0
z
dense_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
░	╨A*
shared_namedense_4/kernel
s
"dense_4/kernel/Read/ReadVariableOpReadVariableOpdense_4/kernel* 
_output_shapes
:
░	╨A*
dtype0
q
dense_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:╨A*
shared_namedense_4/bias
j
 dense_4/bias/Read/ReadVariableOpReadVariableOpdense_4/bias*
_output_shapes	
:╨A*
dtype0
А
conv1d_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:░	а* 
shared_nameconv1d_2/kernel
y
#conv1d_2/kernel/Read/ReadVariableOpReadVariableOpconv1d_2/kernel*$
_output_shapes
:░	а*
dtype0
s
conv1d_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:а*
shared_nameconv1d_2/bias
l
!conv1d_2/bias/Read/ReadVariableOpReadVariableOpconv1d_2/bias*
_output_shapes	
:а*
dtype0
{
dense_6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:аГР*
shared_namedense_6/kernel
t
"dense_6/kernel/Read/ReadVariableOpReadVariableOpdense_6/kernel*!
_output_shapes
:аГР*
dtype0
q
dense_6/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:Р*
shared_namedense_6/bias
j
 dense_6/bias/Read/ReadVariableOpReadVariableOpdense_6/bias*
_output_shapes	
:Р*
dtype0
y
dense_7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	РX*
shared_namedense_7/kernel
r
"dense_7/kernel/Read/ReadVariableOpReadVariableOpdense_7/kernel*
_output_shapes
:	РX*
dtype0
p
dense_7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:X*
shared_namedense_7/bias
i
 dense_7/bias/Read/ReadVariableOpReadVariableOpdense_7/bias*
_output_shapes
:X*
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
Д
conv1d_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:░	а*"
shared_nameconv1d_2/kernel/m
}
%conv1d_2/kernel/m/Read/ReadVariableOpReadVariableOpconv1d_2/kernel/m*$
_output_shapes
:░	а*
dtype0
w
conv1d_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:а* 
shared_nameconv1d_2/bias/m
p
#conv1d_2/bias/m/Read/ReadVariableOpReadVariableOpconv1d_2/bias/m*
_output_shapes	
:а*
dtype0

dense_6/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:аГР*!
shared_namedense_6/kernel/m
x
$dense_6/kernel/m/Read/ReadVariableOpReadVariableOpdense_6/kernel/m*!
_output_shapes
:аГР*
dtype0
u
dense_6/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:Р*
shared_namedense_6/bias/m
n
"dense_6/bias/m/Read/ReadVariableOpReadVariableOpdense_6/bias/m*
_output_shapes	
:Р*
dtype0
}
dense_7/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	РX*!
shared_namedense_7/kernel/m
v
$dense_7/kernel/m/Read/ReadVariableOpReadVariableOpdense_7/kernel/m*
_output_shapes
:	РX*
dtype0
t
dense_7/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:X*
shared_namedense_7/bias/m
m
"dense_7/bias/m/Read/ReadVariableOpReadVariableOpdense_7/bias/m*
_output_shapes
:X*
dtype0
Д
conv1d_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:░	а*"
shared_nameconv1d_2/kernel/v
}
%conv1d_2/kernel/v/Read/ReadVariableOpReadVariableOpconv1d_2/kernel/v*$
_output_shapes
:░	а*
dtype0
w
conv1d_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:а* 
shared_nameconv1d_2/bias/v
p
#conv1d_2/bias/v/Read/ReadVariableOpReadVariableOpconv1d_2/bias/v*
_output_shapes	
:а*
dtype0

dense_6/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:аГР*!
shared_namedense_6/kernel/v
x
$dense_6/kernel/v/Read/ReadVariableOpReadVariableOpdense_6/kernel/v*!
_output_shapes
:аГР*
dtype0
u
dense_6/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:Р*
shared_namedense_6/bias/v
n
"dense_6/bias/v/Read/ReadVariableOpReadVariableOpdense_6/bias/v*
_output_shapes	
:Р*
dtype0
}
dense_7/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	РX*!
shared_namedense_7/kernel/v
v
$dense_7/kernel/v/Read/ReadVariableOpReadVariableOpdense_7/kernel/v*
_output_shapes
:	РX*
dtype0
t
dense_7/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:X*
shared_namedense_7/bias/v
m
"dense_7/bias/v/Read/ReadVariableOpReadVariableOpdense_7/bias/v*
_output_shapes
:X*
dtype0

NoOpNoOp
тM
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*ЭM
valueУMBРM BЙM
╚
layer-0
layer_with_weights-0
layer-1
layer-2
layer_with_weights-1
layer-3
layer_with_weights-2
layer-4
layer_with_weights-3
layer-5
layer_with_weights-4
layer-6
layer_with_weights-5
layer-7
	layer-8

layer-9
layer_with_weights-6
layer-10
layer-11
layer_with_weights-7
layer-12
layer_with_weights-8
layer-13
layer-14
	optimizer

signatures
#_self_saveable_object_factories
trainable_variables
regularization_losses
	variables
	keras_api
%
#_self_saveable_object_factories
Н

kernel
bias
#_self_saveable_object_factories
trainable_variables
regularization_losses
	variables
	keras_api
w
#_self_saveable_object_factories
 trainable_variables
!regularization_losses
"	variables
#	keras_api
Н

$kernel
%bias
#&_self_saveable_object_factories
'trainable_variables
(regularization_losses
)	variables
*	keras_api
Н

+kernel
,bias
#-_self_saveable_object_factories
.trainable_variables
/regularization_losses
0	variables
1	keras_api
Н

2kernel
3bias
#4_self_saveable_object_factories
5trainable_variables
6regularization_losses
7	variables
8	keras_api
Н

9kernel
:bias
#;_self_saveable_object_factories
<trainable_variables
=regularization_losses
>	variables
?	keras_api
Н

@kernel
Abias
#B_self_saveable_object_factories
Ctrainable_variables
Dregularization_losses
E	variables
F	keras_api
w
#G_self_saveable_object_factories
Htrainable_variables
Iregularization_losses
J	variables
K	keras_api
w
#L_self_saveable_object_factories
Mtrainable_variables
Nregularization_losses
O	variables
P	keras_api
Н

Qkernel
Rbias
#S_self_saveable_object_factories
Ttrainable_variables
Uregularization_losses
V	variables
W	keras_api
w
#X_self_saveable_object_factories
Ytrainable_variables
Zregularization_losses
[	variables
\	keras_api
Н

]kernel
^bias
#__self_saveable_object_factories
`trainable_variables
aregularization_losses
b	variables
c	keras_api
Н

dkernel
ebias
#f_self_saveable_object_factories
gtrainable_variables
hregularization_losses
i	variables
j	keras_api
w
#k_self_saveable_object_factories
ltrainable_variables
mregularization_losses
n	variables
o	keras_api
xQm╞Rm╟]m╚^m╔dm╩em╦Qv╠Rv═]v╬^v╧dv╨ev╤
 
 
Ж
0
1
$2
%3
+4
,5
26
37
98
:9
@10
A11
Q12
R13
]14
^15
d16
e17
 
Ж
0
1
$2
%3
+4
,5
26
37
98
:9
@10
A11
Q12
R13
]14
^15
d16
e17
н
player_metrics
trainable_variables
qmetrics
regularization_losses
rnon_trainable_variables

slayers
	variables
tlayer_regularization_losses
 
YW
VARIABLE_VALUEconv1d/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv1d/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE
 

0
1
 

0
1
н
ulayer_metrics
trainable_variables
vmetrics
regularization_losses
wnon_trainable_variables

xlayers
	variables
ylayer_regularization_losses
 
 
 
 
н
zlayer_metrics
 trainable_variables
{metrics
!regularization_losses
|non_trainable_variables

}layers
"	variables
~layer_regularization_losses
XV
VARIABLE_VALUEdense/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
TR
VARIABLE_VALUE
dense/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE
 

$0
%1
 

$0
%1
▒
layer_metrics
'trainable_variables
Аmetrics
(regularization_losses
Бnon_trainable_variables
Вlayers
)	variables
 Гlayer_regularization_losses
ZX
VARIABLE_VALUEdense_1/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_1/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE
 

+0
,1
 

+0
,1
▓
Дlayer_metrics
.trainable_variables
Еmetrics
/regularization_losses
Жnon_trainable_variables
Зlayers
0	variables
 Иlayer_regularization_losses
ZX
VARIABLE_VALUEdense_2/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_2/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE
 

20
31
 

20
31
▓
Йlayer_metrics
5trainable_variables
Кmetrics
6regularization_losses
Лnon_trainable_variables
Мlayers
7	variables
 Нlayer_regularization_losses
ZX
VARIABLE_VALUEdense_3/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_3/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE
 

90
:1
 

90
:1
▓
Оlayer_metrics
<trainable_variables
Пmetrics
=regularization_losses
Рnon_trainable_variables
Сlayers
>	variables
 Тlayer_regularization_losses
ZX
VARIABLE_VALUEdense_4/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_4/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE
 

@0
A1
 

@0
A1
▓
Уlayer_metrics
Ctrainable_variables
Фmetrics
Dregularization_losses
Хnon_trainable_variables
Цlayers
E	variables
 Чlayer_regularization_losses
 
 
 
 
▓
Шlayer_metrics
Htrainable_variables
Щmetrics
Iregularization_losses
Ъnon_trainable_variables
Ыlayers
J	variables
 Ьlayer_regularization_losses
 
 
 
 
▓
Эlayer_metrics
Mtrainable_variables
Юmetrics
Nregularization_losses
Яnon_trainable_variables
аlayers
O	variables
 бlayer_regularization_losses
[Y
VARIABLE_VALUEconv1d_2/kernel6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEconv1d_2/bias4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUE
 

Q0
R1
 

Q0
R1
▓
вlayer_metrics
Ttrainable_variables
гmetrics
Uregularization_losses
дnon_trainable_variables
еlayers
V	variables
 жlayer_regularization_losses
 
 
 
 
▓
зlayer_metrics
Ytrainable_variables
иmetrics
Zregularization_losses
йnon_trainable_variables
кlayers
[	variables
 лlayer_regularization_losses
ZX
VARIABLE_VALUEdense_6/kernel6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_6/bias4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUE
 

]0
^1
 

]0
^1
▓
мlayer_metrics
`trainable_variables
нmetrics
aregularization_losses
оnon_trainable_variables
пlayers
b	variables
 ░layer_regularization_losses
ZX
VARIABLE_VALUEdense_7/kernel6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_7/bias4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUE
 

d0
e1
 

d0
e1
▓
▒layer_metrics
gtrainable_variables
▓metrics
hregularization_losses
│non_trainable_variables
┤layers
i	variables
 ╡layer_regularization_losses
 
 
 
 
▓
╢layer_metrics
ltrainable_variables
╖metrics
mregularization_losses
╕non_trainable_variables
╣layers
n	variables
 ║layer_regularization_losses
 

╗0
╝1
 
n
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
8

╜total

╛count
┐	variables
└	keras_api
I

┴total

┬count
├
_fn_kwargs
─	variables
┼	keras_api
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

╜0
╛1

┐	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

┴0
┬1

─	variables
yw
VARIABLE_VALUEconv1d_2/kernel/mRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEconv1d_2/bias/mPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEdense_6/kernel/mRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEdense_6/bias/mPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEdense_7/kernel/mRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEdense_7/bias/mPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEconv1d_2/kernel/vRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEconv1d_2/bias/vPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEdense_6/kernel/vRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEdense_6/bias/vPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEdense_7/kernel/vRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEdense_7/bias/vPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
Д
serving_default_input_1Placeholder*,
_output_shapes
:         ╨A*
dtype0*!
shape:         ╨A
у
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1conv1d/kernelconv1d/biasdense/kernel
dense/biasdense_1/kerneldense_1/biasdense_2/kerneldense_2/biasdense_3/kerneldense_3/biasdense_4/kerneldense_4/biasconv1d_2/kernelconv1d_2/biasdense_6/kerneldense_6/biasdense_7/kerneldense_7/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2*0,1J 8В *+
f&R$
"__inference_signature_wrapper_2879
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
т
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename!conv1d/kernel/Read/ReadVariableOpconv1d/bias/Read/ReadVariableOp dense/kernel/Read/ReadVariableOpdense/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp dense_1/bias/Read/ReadVariableOp"dense_2/kernel/Read/ReadVariableOp dense_2/bias/Read/ReadVariableOp"dense_3/kernel/Read/ReadVariableOp dense_3/bias/Read/ReadVariableOp"dense_4/kernel/Read/ReadVariableOp dense_4/bias/Read/ReadVariableOp#conv1d_2/kernel/Read/ReadVariableOp!conv1d_2/bias/Read/ReadVariableOp"dense_6/kernel/Read/ReadVariableOp dense_6/bias/Read/ReadVariableOp"dense_7/kernel/Read/ReadVariableOp dense_7/bias/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOp%conv1d_2/kernel/m/Read/ReadVariableOp#conv1d_2/bias/m/Read/ReadVariableOp$dense_6/kernel/m/Read/ReadVariableOp"dense_6/bias/m/Read/ReadVariableOp$dense_7/kernel/m/Read/ReadVariableOp"dense_7/bias/m/Read/ReadVariableOp%conv1d_2/kernel/v/Read/ReadVariableOp#conv1d_2/bias/v/Read/ReadVariableOp$dense_6/kernel/v/Read/ReadVariableOp"dense_6/bias/v/Read/ReadVariableOp$dense_7/kernel/v/Read/ReadVariableOp"dense_7/bias/v/Read/ReadVariableOpConst*/
Tin(
&2$*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *&
f!R
__inference__traced_save_3568
╡
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameconv1d/kernelconv1d/biasdense/kernel
dense/biasdense_1/kerneldense_1/biasdense_2/kerneldense_2/biasdense_3/kerneldense_3/biasdense_4/kerneldense_4/biasconv1d_2/kernelconv1d_2/biasdense_6/kerneldense_6/biasdense_7/kerneldense_7/biastotalcounttotal_1count_1conv1d_2/kernel/mconv1d_2/bias/mdense_6/kernel/mdense_6/bias/mdense_7/kernel/mdense_7/bias/mconv1d_2/kernel/vconv1d_2/bias/vdense_6/kernel/vdense_6/bias/vdense_7/kernel/vdense_7/bias/v*.
Tin'
%2#*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *)
f$R"
 __inference__traced_restore_3680цЦ

ї	
╪
?__inference_dense_layer_call_and_return_conditional_losses_2285

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpР
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         Р2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         Р2

Identity"
identityIdentity:output:0*0
_input_shapes
:         аГ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:Q M
)
_output_shapes
:         аГ
 
_user_specified_nameinputs
┐<
╧
?__inference_model_layer_call_and_return_conditional_losses_2645
input_1
conv1d_2594
conv1d_2596

dense_2600

dense_2602
dense_1_2605
dense_1_2607
dense_2_2610
dense_2_2612
dense_3_2615
dense_3_2617
dense_4_2620
dense_4_2622
conv1d_2_2627
conv1d_2_2629
dense_6_2633
dense_6_2635
dense_7_2638
dense_7_2640
identityИвconv1d/StatefulPartitionedCallв conv1d_2/StatefulPartitionedCallвdense/StatefulPartitionedCallвdense_1/StatefulPartitionedCallвdense_2/StatefulPartitionedCallвdense_3/StatefulPartitionedCallвdense_4/StatefulPartitionedCallвdense_6/StatefulPartitionedCallвdense_7/StatefulPartitionedCallП
conv1d/StatefulPartitionedCallStatefulPartitionedCallinput_1conv1d_2594conv1d_2596*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *I
fDRB
@__inference_conv1d_layer_call_and_return_conditional_losses_22442 
conv1d/StatefulPartitionedCallў
flatten/PartitionedCallPartitionedCall'conv1d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_flatten_layer_call_and_return_conditional_losses_22662
flatten/PartitionedCallЯ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
dense_2600
dense_2602*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_dense_layer_call_and_return_conditional_losses_22852
dense/StatefulPartitionedCallо
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_2605dense_1_2607*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_1_layer_call_and_return_conditional_losses_23122!
dense_1/StatefulPartitionedCall░
dense_2/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0dense_2_2610dense_2_2612*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         <*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_2_layer_call_and_return_conditional_losses_23392!
dense_2/StatefulPartitionedCall▒
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_2615dense_3_2617*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ░	*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_3_layer_call_and_return_conditional_losses_23662!
dense_3/StatefulPartitionedCall▒
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_2620dense_4_2622*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_4_layer_call_and_return_conditional_losses_23932!
dense_4/StatefulPartitionedCallў
dropout/PartitionedCallPartitionedCall(dense_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_24262
dropout/PartitionedCallє
reshape/PartitionedCallPartitionedCall dropout/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_reshape_layer_call_and_return_conditional_losses_24522
reshape/PartitionedCall▓
 conv1d_2/StatefulPartitionedCallStatefulPartitionedCall reshape/PartitionedCall:output:0conv1d_2_2627conv1d_2_2629*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *K
fFRD
B__inference_conv1d_2_layer_call_and_return_conditional_losses_24762"
 conv1d_2/StatefulPartitionedCall∙
dense_5/PartitionedCallPartitionedCall)conv1d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_5_layer_call_and_return_conditional_losses_24982
dense_5/PartitionedCallй
dense_6/StatefulPartitionedCallStatefulPartitionedCall dense_5/PartitionedCall:output:0dense_6_2633dense_6_2635*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_6_layer_call_and_return_conditional_losses_25172!
dense_6/StatefulPartitionedCall░
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_2638dense_7_2640*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_7_layer_call_and_return_conditional_losses_25442!
dense_7/StatefulPartitionedCall№
dropout_2/PartitionedCallPartitionedCall(dense_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *L
fGRE
C__inference_dropout_2_layer_call_and_return_conditional_losses_25772
dropout_2/PartitionedCallж
IdentityIdentity"dropout_2/PartitionedCall:output:0^conv1d/StatefulPartitionedCall!^conv1d_2/StatefulPartitionedCall^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::2@
conv1d/StatefulPartitionedCallconv1d/StatefulPartitionedCall2D
 conv1d_2/StatefulPartitionedCall conv1d_2/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall:U Q
,
_output_shapes
:         ╨A
!
_user_specified_name	input_1
╢?
Ф
?__inference_model_layer_call_and_return_conditional_losses_2702

inputs
conv1d_2651
conv1d_2653

dense_2657

dense_2659
dense_1_2662
dense_1_2664
dense_2_2667
dense_2_2669
dense_3_2672
dense_3_2674
dense_4_2677
dense_4_2679
conv1d_2_2684
conv1d_2_2686
dense_6_2690
dense_6_2692
dense_7_2695
dense_7_2697
identityИвconv1d/StatefulPartitionedCallв conv1d_2/StatefulPartitionedCallвdense/StatefulPartitionedCallвdense_1/StatefulPartitionedCallвdense_2/StatefulPartitionedCallвdense_3/StatefulPartitionedCallвdense_4/StatefulPartitionedCallвdense_6/StatefulPartitionedCallвdense_7/StatefulPartitionedCallвdropout/StatefulPartitionedCallв!dropout_2/StatefulPartitionedCallО
conv1d/StatefulPartitionedCallStatefulPartitionedCallinputsconv1d_2651conv1d_2653*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *I
fDRB
@__inference_conv1d_layer_call_and_return_conditional_losses_22442 
conv1d/StatefulPartitionedCallў
flatten/PartitionedCallPartitionedCall'conv1d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_flatten_layer_call_and_return_conditional_losses_22662
flatten/PartitionedCallЯ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
dense_2657
dense_2659*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_dense_layer_call_and_return_conditional_losses_22852
dense/StatefulPartitionedCallо
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_2662dense_1_2664*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_1_layer_call_and_return_conditional_losses_23122!
dense_1/StatefulPartitionedCall░
dense_2/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0dense_2_2667dense_2_2669*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         <*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_2_layer_call_and_return_conditional_losses_23392!
dense_2/StatefulPartitionedCall▒
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_2672dense_3_2674*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ░	*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_3_layer_call_and_return_conditional_losses_23662!
dense_3/StatefulPartitionedCall▒
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_2677dense_4_2679*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_4_layer_call_and_return_conditional_losses_23932!
dense_4/StatefulPartitionedCallП
dropout/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_24212!
dropout/StatefulPartitionedCall√
reshape/PartitionedCallPartitionedCall(dropout/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_reshape_layer_call_and_return_conditional_losses_24522
reshape/PartitionedCall▓
 conv1d_2/StatefulPartitionedCallStatefulPartitionedCall reshape/PartitionedCall:output:0conv1d_2_2684conv1d_2_2686*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *K
fFRD
B__inference_conv1d_2_layer_call_and_return_conditional_losses_24762"
 conv1d_2/StatefulPartitionedCall∙
dense_5/PartitionedCallPartitionedCall)conv1d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_5_layer_call_and_return_conditional_losses_24982
dense_5/PartitionedCallй
dense_6/StatefulPartitionedCallStatefulPartitionedCall dense_5/PartitionedCall:output:0dense_6_2690dense_6_2692*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_6_layer_call_and_return_conditional_losses_25172!
dense_6/StatefulPartitionedCall░
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_2695dense_7_2697*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_7_layer_call_and_return_conditional_losses_25442!
dense_7/StatefulPartitionedCall╢
!dropout_2/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *L
fGRE
C__inference_dropout_2_layer_call_and_return_conditional_losses_25722#
!dropout_2/StatefulPartitionedCallЇ
IdentityIdentity*dropout_2/StatefulPartitionedCall:output:0^conv1d/StatefulPartitionedCall!^conv1d_2/StatefulPartitionedCall^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dropout/StatefulPartitionedCall"^dropout_2/StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::2@
conv1d/StatefulPartitionedCallconv1d/StatefulPartitionedCall2D
 conv1d_2/StatefulPartitionedCall conv1d_2/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_2/StatefulPartitionedCall!dropout_2/StatefulPartitionedCall:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
╞
a
C__inference_dropout_2_layer_call_and_return_conditional_losses_2577

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:         X2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:         X2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:         X:O K
'
_output_shapes
:         X
 
_user_specified_nameinputs
є

ё
$__inference_model_layer_call_fn_2741
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identityИвStatefulPartitionedCall╠
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_27022
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
,
_output_shapes
:         ╨A
!
_user_specified_name	input_1
Ю
ї
B__inference_conv1d_2_layer_call_and_return_conditional_losses_2476

inputs/
+conv1d_expanddims_1_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpв"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2
conv1d/ExpandDims/dimЧ
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d/ExpandDims║
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dim╣
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d/ExpandDims_1╕
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1dУ
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d/SqueezeН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02
BiasAdd/ReadVariableOpН
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2	
BiasAdd]
ReluReluBiasAdd:output:0*
T0*,
_output_shapes
:         а2
Reluй
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:         а2

Identity"
identityIdentity:output:0*3
_input_shapes"
 :         ╨A::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
█
{
&__inference_dense_2_layer_call_fn_3255

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCallЎ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         <*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_2_layer_call_and_return_conditional_losses_23392
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         <2

Identity"
identityIdentity:output:0*.
_input_shapes
:         ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
▀
]
A__inference_reshape_layer_call_and_return_conditional_losses_2452

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slicee
Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value
B :╨A2
Reshape/shape/1d
Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/2а
Reshape/shapePackstrided_slice:output:0Reshape/shape/1:output:0Reshape/shape/2:output:0*
N*
T0*
_output_shapes
:2
Reshape/shapet
ReshapeReshapeinputsReshape/shape:output:0*
T0*,
_output_shapes
:         ╨A2	
Reshapei
IdentityIdentityReshape:output:0*
T0*,
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
ёЕ
Ї
?__inference_model_layer_call_and_return_conditional_losses_2985

inputs6
2conv1d_conv1d_expanddims_1_readvariableop_resource*
&conv1d_biasadd_readvariableop_resource(
$dense_matmul_readvariableop_resource)
%dense_biasadd_readvariableop_resource*
&dense_1_matmul_readvariableop_resource+
'dense_1_biasadd_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource*
&dense_3_matmul_readvariableop_resource+
'dense_3_biasadd_readvariableop_resource*
&dense_4_matmul_readvariableop_resource+
'dense_4_biasadd_readvariableop_resource8
4conv1d_2_conv1d_expanddims_1_readvariableop_resource,
(conv1d_2_biasadd_readvariableop_resource*
&dense_6_matmul_readvariableop_resource+
'dense_6_biasadd_readvariableop_resource*
&dense_7_matmul_readvariableop_resource+
'dense_7_biasadd_readvariableop_resource
identityИвconv1d/BiasAdd/ReadVariableOpв)conv1d/conv1d/ExpandDims_1/ReadVariableOpвconv1d_2/BiasAdd/ReadVariableOpв+conv1d_2/conv1d/ExpandDims_1/ReadVariableOpвdense/BiasAdd/ReadVariableOpвdense/MatMul/ReadVariableOpвdense_1/BiasAdd/ReadVariableOpвdense_1/MatMul/ReadVariableOpвdense_2/BiasAdd/ReadVariableOpвdense_2/MatMul/ReadVariableOpвdense_3/BiasAdd/ReadVariableOpвdense_3/MatMul/ReadVariableOpвdense_4/BiasAdd/ReadVariableOpвdense_4/MatMul/ReadVariableOpвdense_6/BiasAdd/ReadVariableOpвdense_6/MatMul/ReadVariableOpвdense_7/BiasAdd/ReadVariableOpвdense_7/MatMul/ReadVariableOpЗ
conv1d/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2
conv1d/conv1d/ExpandDims/dimм
conv1d/conv1d/ExpandDims
ExpandDimsinputs%conv1d/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d/conv1d/ExpandDims╧
)conv1d/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv1d_conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02+
)conv1d/conv1d/ExpandDims_1/ReadVariableOpВ
conv1d/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv1d/conv1d/ExpandDims_1/dim╒
conv1d/conv1d/ExpandDims_1
ExpandDims1conv1d/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv1d/conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d/conv1d/ExpandDims_1╘
conv1d/conv1dConv2D!conv1d/conv1d/ExpandDims:output:0#conv1d/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1d/conv1dи
conv1d/conv1d/SqueezeSqueezeconv1d/conv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d/conv1d/Squeezeв
conv1d/BiasAdd/ReadVariableOpReadVariableOp&conv1d_biasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02
conv1d/BiasAdd/ReadVariableOpй
conv1d/BiasAddBiasAddconv1d/conv1d/Squeeze:output:0%conv1d/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2
conv1d/BiasAddr
conv1d/ReluReluconv1d/BiasAdd:output:0*
T0*,
_output_shapes
:         а2
conv1d/Reluo
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
flatten/ConstФ
flatten/ReshapeReshapeconv1d/Relu:activations:0flatten/Const:output:0*
T0*)
_output_shapes
:         аГ2
flatten/Reshapeв
dense/MatMul/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
dense/MatMul/ReadVariableOpШ
dense/MatMulMatMulflatten/Reshape:output:0#dense/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense/MatMulЯ
dense/BiasAdd/ReadVariableOpReadVariableOp%dense_biasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02
dense/BiasAdd/ReadVariableOpЪ
dense/BiasAddBiasAdddense/MatMul:product:0$dense/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense/BiasAddk

dense/ReluReludense/BiasAdd:output:0*
T0*(
_output_shapes
:         Р2

dense/Reluж
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	Р*
dtype02
dense_1/MatMul/ReadVariableOpЭ
dense_1/MatMulMatMuldense/Relu:activations:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
dense_1/MatMulд
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_1/BiasAdd/ReadVariableOpб
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
dense_1/BiasAddp
dense_1/TanhTanhdense_1/BiasAdd:output:0*
T0*'
_output_shapes
:         2
dense_1/Tanhе
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:<*
dtype02
dense_2/MatMul/ReadVariableOpХ
dense_2/MatMulMatMuldense_1/Tanh:y:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
dense_2/MatMulд
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:<*
dtype02 
dense_2/BiasAdd/ReadVariableOpб
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
dense_2/BiasAddp
dense_2/TanhTanhdense_2/BiasAdd:output:0*
T0*'
_output_shapes
:         <2
dense_2/Tanhж
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes
:	<░	*
dtype02
dense_3/MatMul/ReadVariableOpЦ
dense_3/MatMulMatMuldense_2/Tanh:y:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
dense_3/MatMulе
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes	
:░	*
dtype02 
dense_3/BiasAdd/ReadVariableOpв
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
dense_3/BiasAddq
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*(
_output_shapes
:         ░	2
dense_3/Reluз
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource* 
_output_shapes
:
░	╨A*
dtype02
dense_4/MatMul/ReadVariableOpа
dense_4/MatMulMatMuldense_3/Relu:activations:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
dense_4/MatMulе
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes	
:╨A*
dtype02 
dense_4/BiasAdd/ReadVariableOpв
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
dense_4/BiasAddq
dense_4/ReluReludense_4/BiasAdd:output:0*
T0*(
_output_shapes
:         ╨A2
dense_4/Relus
dropout/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
dropout/dropout/Constа
dropout/dropout/MulMuldense_4/Relu:activations:0dropout/dropout/Const:output:0*
T0*(
_output_shapes
:         ╨A2
dropout/dropout/Mulx
dropout/dropout/ShapeShapedense_4/Relu:activations:0*
T0*
_output_shapes
:2
dropout/dropout/Shape═
,dropout/dropout/random_uniform/RandomUniformRandomUniformdropout/dropout/Shape:output:0*
T0*(
_output_shapes
:         ╨A*
dtype02.
,dropout/dropout/random_uniform/RandomUniformЕ
dropout/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dropout/dropout/GreaterEqual/y▀
dropout/dropout/GreaterEqualGreaterEqual5dropout/dropout/random_uniform/RandomUniform:output:0'dropout/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:         ╨A2
dropout/dropout/GreaterEqualШ
dropout/dropout/CastCast dropout/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:         ╨A2
dropout/dropout/CastЫ
dropout/dropout/Mul_1Muldropout/dropout/Mul:z:0dropout/dropout/Cast:y:0*
T0*(
_output_shapes
:         ╨A2
dropout/dropout/Mul_1g
reshape/ShapeShapedropout/dropout/Mul_1:z:0*
T0*
_output_shapes
:2
reshape/ShapeД
reshape/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
reshape/strided_slice/stackИ
reshape/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
reshape/strided_slice/stack_1И
reshape/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
reshape/strided_slice/stack_2Т
reshape/strided_sliceStridedSlicereshape/Shape:output:0$reshape/strided_slice/stack:output:0&reshape/strided_slice/stack_1:output:0&reshape/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
reshape/strided_sliceu
reshape/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value
B :╨A2
reshape/Reshape/shape/1t
reshape/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
reshape/Reshape/shape/2╚
reshape/Reshape/shapePackreshape/strided_slice:output:0 reshape/Reshape/shape/1:output:0 reshape/Reshape/shape/2:output:0*
N*
T0*
_output_shapes
:2
reshape/Reshape/shapeЯ
reshape/ReshapeReshapedropout/dropout/Mul_1:z:0reshape/Reshape/shape:output:0*
T0*,
_output_shapes
:         ╨A2
reshape/ReshapeЛ
conv1d_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2 
conv1d_2/conv1d/ExpandDims/dim─
conv1d_2/conv1d/ExpandDims
ExpandDimsreshape/Reshape:output:0'conv1d_2/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d_2/conv1d/ExpandDims╒
+conv1d_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp4conv1d_2_conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02-
+conv1d_2/conv1d/ExpandDims_1/ReadVariableOpЖ
 conv1d_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2"
 conv1d_2/conv1d/ExpandDims_1/dim▌
conv1d_2/conv1d/ExpandDims_1
ExpandDims3conv1d_2/conv1d/ExpandDims_1/ReadVariableOp:value:0)conv1d_2/conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d_2/conv1d/ExpandDims_1▄
conv1d_2/conv1dConv2D#conv1d_2/conv1d/ExpandDims:output:0%conv1d_2/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1d_2/conv1dо
conv1d_2/conv1d/SqueezeSqueezeconv1d_2/conv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d_2/conv1d/Squeezeи
conv1d_2/BiasAdd/ReadVariableOpReadVariableOp(conv1d_2_biasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02!
conv1d_2/BiasAdd/ReadVariableOp▒
conv1d_2/BiasAddBiasAdd conv1d_2/conv1d/Squeeze:output:0'conv1d_2/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2
conv1d_2/BiasAddx
conv1d_2/ReluReluconv1d_2/BiasAdd:output:0*
T0*,
_output_shapes
:         а2
conv1d_2/Reluo
dense_5/ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
dense_5/ConstЦ
dense_5/ReshapeReshapeconv1d_2/Relu:activations:0dense_5/Const:output:0*
T0*)
_output_shapes
:         аГ2
dense_5/Reshapeи
dense_6/MatMul/ReadVariableOpReadVariableOp&dense_6_matmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
dense_6/MatMul/ReadVariableOpЮ
dense_6/MatMulMatMuldense_5/Reshape:output:0%dense_6/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense_6/MatMulе
dense_6/BiasAdd/ReadVariableOpReadVariableOp'dense_6_biasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02 
dense_6/BiasAdd/ReadVariableOpв
dense_6/BiasAddBiasAdddense_6/MatMul:product:0&dense_6/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense_6/BiasAddq
dense_6/ReluReludense_6/BiasAdd:output:0*
T0*(
_output_shapes
:         Р2
dense_6/Reluж
dense_7/MatMul/ReadVariableOpReadVariableOp&dense_7_matmul_readvariableop_resource*
_output_shapes
:	РX*
dtype02
dense_7/MatMul/ReadVariableOpЯ
dense_7/MatMulMatMuldense_6/Relu:activations:0%dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
dense_7/MatMulд
dense_7/BiasAdd/ReadVariableOpReadVariableOp'dense_7_biasadd_readvariableop_resource*
_output_shapes
:X*
dtype02 
dense_7/BiasAdd/ReadVariableOpб
dense_7/BiasAddBiasAdddense_7/MatMul:product:0&dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
dense_7/BiasAddy
dense_7/SoftmaxSoftmaxdense_7/BiasAdd:output:0*
T0*'
_output_shapes
:         X2
dense_7/Softmaxw
dropout_2/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
dropout_2/dropout/Constд
dropout_2/dropout/MulMuldense_7/Softmax:softmax:0 dropout_2/dropout/Const:output:0*
T0*'
_output_shapes
:         X2
dropout_2/dropout/Mul{
dropout_2/dropout/ShapeShapedense_7/Softmax:softmax:0*
T0*
_output_shapes
:2
dropout_2/dropout/Shape╥
.dropout_2/dropout/random_uniform/RandomUniformRandomUniform dropout_2/dropout/Shape:output:0*
T0*'
_output_shapes
:         X*
dtype020
.dropout_2/dropout/random_uniform/RandomUniformЙ
 dropout_2/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dropout_2/dropout/GreaterEqual/yц
dropout_2/dropout/GreaterEqualGreaterEqual7dropout_2/dropout/random_uniform/RandomUniform:output:0)dropout_2/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:         X2 
dropout_2/dropout/GreaterEqualЭ
dropout_2/dropout/CastCast"dropout_2/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:         X2
dropout_2/dropout/Castв
dropout_2/dropout/Mul_1Muldropout_2/dropout/Mul:z:0dropout_2/dropout/Cast:y:0*
T0*'
_output_shapes
:         X2
dropout_2/dropout/Mul_1╬
IdentityIdentitydropout_2/dropout/Mul_1:z:0^conv1d/BiasAdd/ReadVariableOp*^conv1d/conv1d/ExpandDims_1/ReadVariableOp ^conv1d_2/BiasAdd/ReadVariableOp,^conv1d_2/conv1d/ExpandDims_1/ReadVariableOp^dense/BiasAdd/ReadVariableOp^dense/MatMul/ReadVariableOp^dense_1/BiasAdd/ReadVariableOp^dense_1/MatMul/ReadVariableOp^dense_2/BiasAdd/ReadVariableOp^dense_2/MatMul/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp^dense_6/BiasAdd/ReadVariableOp^dense_6/MatMul/ReadVariableOp^dense_7/BiasAdd/ReadVariableOp^dense_7/MatMul/ReadVariableOp*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::2>
conv1d/BiasAdd/ReadVariableOpconv1d/BiasAdd/ReadVariableOp2V
)conv1d/conv1d/ExpandDims_1/ReadVariableOp)conv1d/conv1d/ExpandDims_1/ReadVariableOp2B
conv1d_2/BiasAdd/ReadVariableOpconv1d_2/BiasAdd/ReadVariableOp2Z
+conv1d_2/conv1d/ExpandDims_1/ReadVariableOp+conv1d_2/conv1d/ExpandDims_1/ReadVariableOp2<
dense/BiasAdd/ReadVariableOpdense/BiasAdd/ReadVariableOp2:
dense/MatMul/ReadVariableOpdense/MatMul/ReadVariableOp2@
dense_1/BiasAdd/ReadVariableOpdense_1/BiasAdd/ReadVariableOp2>
dense_1/MatMul/ReadVariableOpdense_1/MatMul/ReadVariableOp2@
dense_2/BiasAdd/ReadVariableOpdense_2/BiasAdd/ReadVariableOp2>
dense_2/MatMul/ReadVariableOpdense_2/MatMul/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp2@
dense_6/BiasAdd/ReadVariableOpdense_6/BiasAdd/ReadVariableOp2>
dense_6/MatMul/ReadVariableOpdense_6/MatMul/ReadVariableOp2@
dense_7/BiasAdd/ReadVariableOpdense_7/BiasAdd/ReadVariableOp2>
dense_7/MatMul/ReadVariableOpdense_7/MatMul/ReadVariableOp:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
 

b
C__inference_dropout_2_layer_call_and_return_conditional_losses_2572

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:         X2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape┤
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:         X*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *    2
dropout/GreaterEqual/y╛
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:         X2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:         X2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:         X2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*&
_input_shapes
:         X:O K
'
_output_shapes
:         X
 
_user_specified_nameinputs
с
{
&__inference_dense_6_layer_call_fn_3396

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCallў
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_6_layer_call_and_return_conditional_losses_25172
StatefulPartitionedCallП
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:         Р2

Identity"
identityIdentity:output:0*0
_input_shapes
:         аГ::22
StatefulPartitionedCallStatefulPartitionedCall:Q M
)
_output_shapes
:         аГ
 
_user_specified_nameinputs
Ё

Ё
$__inference_model_layer_call_fn_3118

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identityИвStatefulPartitionedCall╦
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_27022
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
ё	
┌
A__inference_dense_3_layer_call_and_return_conditional_losses_2366

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	<░	*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:░	*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         ░	2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         ░	2

Identity"
identityIdentity:output:0*.
_input_shapes
:         <::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         <
 
_user_specified_nameinputs
ў	
┌
A__inference_dense_6_layer_call_and_return_conditional_losses_3387

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpР
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         Р2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         Р2

Identity"
identityIdentity:output:0*0
_input_shapes
:         аГ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:Q M
)
_output_shapes
:         аГ
 
_user_specified_nameinputs
Ў	
┌
A__inference_dense_7_layer_call_and_return_conditional_losses_3407

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	РX*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:X*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2	
BiasAdda
SoftmaxSoftmaxBiasAdd:output:0*
T0*'
_output_shapes
:         X2	
SoftmaxЦ
IdentityIdentitySoftmax:softmax:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*/
_input_shapes
:         Р::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:         Р
 
_user_specified_nameinputs
╡
]
A__inference_dense_5_layer_call_and_return_conditional_losses_3371

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
Consti
ReshapeReshapeinputsConst:output:0*
T0*)
_output_shapes
:         аГ2	
Reshapef
IdentityIdentityReshape:output:0*
T0*)
_output_shapes
:         аГ2

Identity"
identityIdentity:output:0*+
_input_shapes
:         а:T P
,
_output_shapes
:         а
 
_user_specified_nameinputs
ї	
╪
?__inference_dense_layer_call_and_return_conditional_losses_3206

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpР
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         Р2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         Р2

Identity"
identityIdentity:output:0*0
_input_shapes
:         аГ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:Q M
)
_output_shapes
:         аГ
 
_user_specified_nameinputs
ў	
┌
A__inference_dense_6_layer_call_and_return_conditional_losses_2517

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpР
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         Р2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         Р2

Identity"
identityIdentity:output:0*0
_input_shapes
:         аГ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:Q M
)
_output_shapes
:         аГ
 
_user_specified_nameinputs
Я
B
&__inference_flatten_layer_call_fn_3195

inputs
identity╞
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_flatten_layer_call_and_return_conditional_losses_22662
PartitionedCalln
IdentityIdentityPartitionedCall:output:0*
T0*)
_output_shapes
:         аГ2

Identity"
identityIdentity:output:0*+
_input_shapes
:         а:T P
,
_output_shapes
:         а
 
_user_specified_nameinputs
Х
B
&__inference_dropout_layer_call_fn_3322

inputs
identity┼
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_24262
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
ФH
В
__inference__traced_save_3568
file_prefix,
(savev2_conv1d_kernel_read_readvariableop*
&savev2_conv1d_bias_read_readvariableop+
'savev2_dense_kernel_read_readvariableop)
%savev2_dense_bias_read_readvariableop-
)savev2_dense_1_kernel_read_readvariableop+
'savev2_dense_1_bias_read_readvariableop-
)savev2_dense_2_kernel_read_readvariableop+
'savev2_dense_2_bias_read_readvariableop-
)savev2_dense_3_kernel_read_readvariableop+
'savev2_dense_3_bias_read_readvariableop-
)savev2_dense_4_kernel_read_readvariableop+
'savev2_dense_4_bias_read_readvariableop.
*savev2_conv1d_2_kernel_read_readvariableop,
(savev2_conv1d_2_bias_read_readvariableop-
)savev2_dense_6_kernel_read_readvariableop+
'savev2_dense_6_bias_read_readvariableop-
)savev2_dense_7_kernel_read_readvariableop+
'savev2_dense_7_bias_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop0
,savev2_conv1d_2_kernel_m_read_readvariableop.
*savev2_conv1d_2_bias_m_read_readvariableop/
+savev2_dense_6_kernel_m_read_readvariableop-
)savev2_dense_6_bias_m_read_readvariableop/
+savev2_dense_7_kernel_m_read_readvariableop-
)savev2_dense_7_bias_m_read_readvariableop0
,savev2_conv1d_2_kernel_v_read_readvariableop.
*savev2_conv1d_2_bias_v_read_readvariableop/
+savev2_dense_6_kernel_v_read_readvariableop-
)savev2_dense_6_bias_v_read_readvariableop/
+savev2_dense_7_kernel_v_read_readvariableop-
)savev2_dense_7_bias_v_read_readvariableop
savev2_const

identity_1ИвMergeV2CheckpointsП
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Constl
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part2	
Const_1Л
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shardж
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename╗
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:#*
dtype0*═
value├B└#B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_names╬
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:#*
dtype0*Y
valuePBN#B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slicesь
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0(savev2_conv1d_kernel_read_readvariableop&savev2_conv1d_bias_read_readvariableop'savev2_dense_kernel_read_readvariableop%savev2_dense_bias_read_readvariableop)savev2_dense_1_kernel_read_readvariableop'savev2_dense_1_bias_read_readvariableop)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop)savev2_dense_3_kernel_read_readvariableop'savev2_dense_3_bias_read_readvariableop)savev2_dense_4_kernel_read_readvariableop'savev2_dense_4_bias_read_readvariableop*savev2_conv1d_2_kernel_read_readvariableop(savev2_conv1d_2_bias_read_readvariableop)savev2_dense_6_kernel_read_readvariableop'savev2_dense_6_bias_read_readvariableop)savev2_dense_7_kernel_read_readvariableop'savev2_dense_7_bias_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop,savev2_conv1d_2_kernel_m_read_readvariableop*savev2_conv1d_2_bias_m_read_readvariableop+savev2_dense_6_kernel_m_read_readvariableop)savev2_dense_6_bias_m_read_readvariableop+savev2_dense_7_kernel_m_read_readvariableop)savev2_dense_7_bias_m_read_readvariableop,savev2_conv1d_2_kernel_v_read_readvariableop*savev2_conv1d_2_bias_v_read_readvariableop+savev2_dense_6_kernel_v_read_readvariableop)savev2_dense_6_bias_v_read_readvariableop+savev2_dense_7_kernel_v_read_readvariableop)savev2_dense_7_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *1
dtypes'
%2#2
SaveV2║
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixesб
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identitym

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*╞
_input_shapes┤
▒: :░	а:а:аГР:Р:	Р::<:<:	<░	:░	:
░	╨A:╨A:░	а:а:аГР:Р:	РX:X: : : : :░	а:а:аГР:Р:	РX:X:░	а:а:аГР:Р:	РX:X: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:*&
$
_output_shapes
:░	а:!

_output_shapes	
:а:'#
!
_output_shapes
:аГР:!

_output_shapes	
:Р:%!

_output_shapes
:	Р: 

_output_shapes
::$ 

_output_shapes

:<: 

_output_shapes
:<:%	!

_output_shapes
:	<░	:!


_output_shapes	
:░	:&"
 
_output_shapes
:
░	╨A:!

_output_shapes	
:╨A:*&
$
_output_shapes
:░	а:!

_output_shapes	
:а:'#
!
_output_shapes
:аГР:!

_output_shapes	
:Р:%!

_output_shapes
:	РX: 

_output_shapes
:X:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :*&
$
_output_shapes
:░	а:!

_output_shapes	
:а:'#
!
_output_shapes
:аГР:!

_output_shapes	
:Р:%!

_output_shapes
:	РX: 

_output_shapes
:X:*&
$
_output_shapes
:░	а:!

_output_shapes	
:а:'#
!
_output_shapes
:аГР:! 

_output_shapes	
:Р:%!!

_output_shapes
:	РX: "

_output_shapes
:X:#

_output_shapes
: 
Х
D
(__inference_dropout_2_layer_call_fn_3443

inputs
identity╞
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *L
fGRE
C__inference_dropout_2_layer_call_and_return_conditional_losses_25772
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*&
_input_shapes
:         X:O K
'
_output_shapes
:         X
 
_user_specified_nameinputs
Ё

Ё
$__inference_model_layer_call_fn_3159

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identityИвStatefulPartitionedCall╦
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_27972
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
╞
a
C__inference_dropout_2_layer_call_and_return_conditional_losses_3433

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:         X2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:         X2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:         X:O K
'
_output_shapes
:         X
 
_user_specified_nameinputs
▀
]
A__inference_reshape_layer_call_and_return_conditional_losses_3335

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slicee
Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value
B :╨A2
Reshape/shape/1d
Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/2а
Reshape/shapePackstrided_slice:output:0Reshape/shape/1:output:0Reshape/shape/2:output:0*
N*
T0*
_output_shapes
:2
Reshape/shapet
ReshapeReshapeinputsReshape/shape:output:0*
T0*,
_output_shapes
:         ╨A2	
Reshapei
IdentityIdentityReshape:output:0*
T0*,
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
ф	
┌
A__inference_dense_1_layer_call_and_return_conditional_losses_3226

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	Р*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:         2
TanhН
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:         2

Identity"
identityIdentity:output:0*/
_input_shapes
:         Р::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:         Р
 
_user_specified_nameinputs
э
z
%__inference_conv1d_layer_call_fn_3184

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCall·
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *I
fDRB
@__inference_conv1d_layer_call_and_return_conditional_losses_22442
StatefulPartitionedCallУ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*,
_output_shapes
:         а2

Identity"
identityIdentity:output:0*3
_input_shapes"
 :         ╨A::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
Ж
`
A__inference_dropout_layer_call_and_return_conditional_losses_2421

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:         ╨A2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape╡
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:         ╨A*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *    2
dropout/GreaterEqual/y┐
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:         ╨A2
dropout/GreaterEqualА
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:         ╨A2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:         ╨A2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
Ь
є
@__inference_conv1d_layer_call_and_return_conditional_losses_2244

inputs/
+conv1d_expanddims_1_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpв"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2
conv1d/ExpandDims/dimЧ
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d/ExpandDims║
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dim╣
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d/ExpandDims_1╕
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1dУ
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d/SqueezeН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02
BiasAdd/ReadVariableOpН
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2	
BiasAdd]
ReluReluBiasAdd:output:0*
T0*,
_output_shapes
:         а2
Reluй
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:         а2

Identity"
identityIdentity:output:0*3
_input_shapes"
 :         ╨A::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
╚
_
A__inference_dropout_layer_call_and_return_conditional_losses_3312

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:         ╨A2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:         ╨A2

Identity_1"!

identity_1Identity_1:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
╝<
╬
?__inference_model_layer_call_and_return_conditional_losses_2797

inputs
conv1d_2746
conv1d_2748

dense_2752

dense_2754
dense_1_2757
dense_1_2759
dense_2_2762
dense_2_2764
dense_3_2767
dense_3_2769
dense_4_2772
dense_4_2774
conv1d_2_2779
conv1d_2_2781
dense_6_2785
dense_6_2787
dense_7_2790
dense_7_2792
identityИвconv1d/StatefulPartitionedCallв conv1d_2/StatefulPartitionedCallвdense/StatefulPartitionedCallвdense_1/StatefulPartitionedCallвdense_2/StatefulPartitionedCallвdense_3/StatefulPartitionedCallвdense_4/StatefulPartitionedCallвdense_6/StatefulPartitionedCallвdense_7/StatefulPartitionedCallО
conv1d/StatefulPartitionedCallStatefulPartitionedCallinputsconv1d_2746conv1d_2748*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *I
fDRB
@__inference_conv1d_layer_call_and_return_conditional_losses_22442 
conv1d/StatefulPartitionedCallў
flatten/PartitionedCallPartitionedCall'conv1d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_flatten_layer_call_and_return_conditional_losses_22662
flatten/PartitionedCallЯ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
dense_2752
dense_2754*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_dense_layer_call_and_return_conditional_losses_22852
dense/StatefulPartitionedCallо
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_2757dense_1_2759*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_1_layer_call_and_return_conditional_losses_23122!
dense_1/StatefulPartitionedCall░
dense_2/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0dense_2_2762dense_2_2764*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         <*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_2_layer_call_and_return_conditional_losses_23392!
dense_2/StatefulPartitionedCall▒
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_2767dense_3_2769*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ░	*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_3_layer_call_and_return_conditional_losses_23662!
dense_3/StatefulPartitionedCall▒
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_2772dense_4_2774*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_4_layer_call_and_return_conditional_losses_23932!
dense_4/StatefulPartitionedCallў
dropout/PartitionedCallPartitionedCall(dense_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_24262
dropout/PartitionedCallє
reshape/PartitionedCallPartitionedCall dropout/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_reshape_layer_call_and_return_conditional_losses_24522
reshape/PartitionedCall▓
 conv1d_2/StatefulPartitionedCallStatefulPartitionedCall reshape/PartitionedCall:output:0conv1d_2_2779conv1d_2_2781*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *K
fFRD
B__inference_conv1d_2_layer_call_and_return_conditional_losses_24762"
 conv1d_2/StatefulPartitionedCall∙
dense_5/PartitionedCallPartitionedCall)conv1d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_5_layer_call_and_return_conditional_losses_24982
dense_5/PartitionedCallй
dense_6/StatefulPartitionedCallStatefulPartitionedCall dense_5/PartitionedCall:output:0dense_6_2785dense_6_2787*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_6_layer_call_and_return_conditional_losses_25172!
dense_6/StatefulPartitionedCall░
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_2790dense_7_2792*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_7_layer_call_and_return_conditional_losses_25442!
dense_7/StatefulPartitionedCall№
dropout_2/PartitionedCallPartitionedCall(dense_7/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *L
fGRE
C__inference_dropout_2_layer_call_and_return_conditional_losses_25772
dropout_2/PartitionedCallж
IdentityIdentity"dropout_2/PartitionedCall:output:0^conv1d/StatefulPartitionedCall!^conv1d_2/StatefulPartitionedCall^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::2@
conv1d/StatefulPartitionedCallconv1d/StatefulPartitionedCall2D
 conv1d_2/StatefulPartitionedCall conv1d_2/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
ф	
┌
A__inference_dense_1_layer_call_and_return_conditional_losses_2312

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	Р*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:         2
TanhН
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:         2

Identity"
identityIdentity:output:0*/
_input_shapes
:         Р::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:         Р
 
_user_specified_nameinputs
╚
_
A__inference_dropout_layer_call_and_return_conditional_losses_2426

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:         ╨A2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:         ╨A2

Identity_1"!

identity_1Identity_1:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
Я
B
&__inference_dense_5_layer_call_fn_3376

inputs
identity╞
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_5_layer_call_and_return_conditional_losses_24982
PartitionedCalln
IdentityIdentityPartitionedCall:output:0*
T0*)
_output_shapes
:         аГ2

Identity"
identityIdentity:output:0*+
_input_shapes
:         а:T P
,
_output_shapes
:         а
 
_user_specified_nameinputs
с	
┌
A__inference_dense_2_layer_call_and_return_conditional_losses_3246

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpН
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:<*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:<*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:         <2
TanhН
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:         <2

Identity"
identityIdentity:output:0*.
_input_shapes
:         ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
╡
]
A__inference_dense_5_layer_call_and_return_conditional_losses_2498

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
Consti
ReshapeReshapeinputsConst:output:0*
T0*)
_output_shapes
:         аГ2	
Reshapef
IdentityIdentityReshape:output:0*
T0*)
_output_shapes
:         аГ2

Identity"
identityIdentity:output:0*+
_input_shapes
:         а:T P
,
_output_shapes
:         а
 
_user_specified_nameinputs
Ю
ї
B__inference_conv1d_2_layer_call_and_return_conditional_losses_3356

inputs/
+conv1d_expanddims_1_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpв"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2
conv1d/ExpandDims/dimЧ
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d/ExpandDims║
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dim╣
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d/ExpandDims_1╕
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1dУ
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d/SqueezeН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02
BiasAdd/ReadVariableOpН
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2	
BiasAdd]
ReluReluBiasAdd:output:0*
T0*,
_output_shapes
:         а2
Reluй
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:         а2

Identity"
identityIdentity:output:0*3
_input_shapes"
 :         ╨A::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
с	
┌
A__inference_dense_2_layer_call_and_return_conditional_losses_2339

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpН
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:<*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:<*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2	
BiasAddX
TanhTanhBiasAdd:output:0*
T0*'
_output_shapes
:         <2
TanhН
IdentityIdentityTanh:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:         <2

Identity"
identityIdentity:output:0*.
_input_shapes
:         ::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         
 
_user_specified_nameinputs
 

b
C__inference_dropout_2_layer_call_and_return_conditional_losses_3428

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:         X2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape┤
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:         X*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *    2
dropout/GreaterEqual/y╛
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:         X2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:         X2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:         X2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*&
_input_shapes
:         X:O K
'
_output_shapes
:         X
 
_user_specified_nameinputs
▌
{
&__inference_dense_3_layer_call_fn_3275

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCallў
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ░	*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_3_layer_call_and_return_conditional_losses_23662
StatefulPartitionedCallП
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:         ░	2

Identity"
identityIdentity:output:0*.
_input_shapes
:         <::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         <
 
_user_specified_nameinputs
╤

я
"__inference_signature_wrapper_2879
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identityИвStatefulPartitionedCallм
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2*0,1J 8В *(
f#R!
__inference__wrapped_model_22242
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
,
_output_shapes
:         ╨A
!
_user_specified_name	input_1
Ь
є
@__inference_conv1d_layer_call_and_return_conditional_losses_3175

inputs/
+conv1d_expanddims_1_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpв"conv1d/ExpandDims_1/ReadVariableOpy
conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2
conv1d/ExpandDims/dimЧ
conv1d/ExpandDims
ExpandDimsinputsconv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d/ExpandDims║
"conv1d/ExpandDims_1/ReadVariableOpReadVariableOp+conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02$
"conv1d/ExpandDims_1/ReadVariableOpt
conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2
conv1d/ExpandDims_1/dim╣
conv1d/ExpandDims_1
ExpandDims*conv1d/ExpandDims_1/ReadVariableOp:value:0 conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d/ExpandDims_1╕
conv1dConv2Dconv1d/ExpandDims:output:0conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1dУ
conv1d/SqueezeSqueezeconv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d/SqueezeН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02
BiasAdd/ReadVariableOpН
BiasAddBiasAddconv1d/Squeeze:output:0BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2	
BiasAdd]
ReluReluBiasAdd:output:0*
T0*,
_output_shapes
:         а2
Reluй
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp#^conv1d/ExpandDims_1/ReadVariableOp*
T0*,
_output_shapes
:         а2

Identity"
identityIdentity:output:0*3
_input_shapes"
 :         ╨A::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2H
"conv1d/ExpandDims_1/ReadVariableOp"conv1d/ExpandDims_1/ReadVariableOp:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
є

ё
$__inference_model_layer_call_fn_2836
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16
identityИвStatefulPartitionedCall╠
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*4
_read_only_resource_inputs
	
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_27972
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
,
_output_shapes
:         ╨A
!
_user_specified_name	input_1
╡
]
A__inference_flatten_layer_call_and_return_conditional_losses_3190

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
Consti
ReshapeReshapeinputsConst:output:0*
T0*)
_output_shapes
:         аГ2	
Reshapef
IdentityIdentityReshape:output:0*
T0*)
_output_shapes
:         аГ2

Identity"
identityIdentity:output:0*+
_input_shapes
:         а:T P
,
_output_shapes
:         а
 
_user_specified_nameinputs
Ў	
┌
A__inference_dense_7_layer_call_and_return_conditional_losses_2544

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	РX*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:X*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2	
BiasAdda
SoftmaxSoftmaxBiasAdd:output:0*
T0*'
_output_shapes
:         X2	
SoftmaxЦ
IdentityIdentitySoftmax:softmax:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*/
_input_shapes
:         Р::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:         Р
 
_user_specified_nameinputs
Ї	
┌
A__inference_dense_4_layer_call_and_return_conditional_losses_3286

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpП
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
░	╨A*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:╨A*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         ╨A2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*/
_input_shapes
:         ░	::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:         ░	
 
_user_specified_nameinputs
Ї	
┌
A__inference_dense_4_layer_call_and_return_conditional_losses_2393

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpП
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
░	╨A*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:╨A*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         ╨A2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*/
_input_shapes
:         ░	::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:         ░	
 
_user_specified_nameinputs
б
a
(__inference_dropout_2_layer_call_fn_3438

inputs
identityИвStatefulPartitionedCall▐
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *L
fGRE
C__inference_dropout_2_layer_call_and_return_conditional_losses_25722
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*&
_input_shapes
:         X22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         X
 
_user_specified_nameinputs
▌
{
&__inference_dense_7_layer_call_fn_3416

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCallЎ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_7_layer_call_and_return_conditional_losses_25442
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*/
_input_shapes
:         Р::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         Р
 
_user_specified_nameinputs
ё
|
'__inference_conv1d_2_layer_call_fn_3365

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCall№
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *K
fFRD
B__inference_conv1d_2_layer_call_and_return_conditional_losses_24762
StatefulPartitionedCallУ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*,
_output_shapes
:         а2

Identity"
identityIdentity:output:0*3
_input_shapes"
 :         ╨A::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
┼s
Ї
?__inference_model_layer_call_and_return_conditional_losses_3077

inputs6
2conv1d_conv1d_expanddims_1_readvariableop_resource*
&conv1d_biasadd_readvariableop_resource(
$dense_matmul_readvariableop_resource)
%dense_biasadd_readvariableop_resource*
&dense_1_matmul_readvariableop_resource+
'dense_1_biasadd_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource*
&dense_3_matmul_readvariableop_resource+
'dense_3_biasadd_readvariableop_resource*
&dense_4_matmul_readvariableop_resource+
'dense_4_biasadd_readvariableop_resource8
4conv1d_2_conv1d_expanddims_1_readvariableop_resource,
(conv1d_2_biasadd_readvariableop_resource*
&dense_6_matmul_readvariableop_resource+
'dense_6_biasadd_readvariableop_resource*
&dense_7_matmul_readvariableop_resource+
'dense_7_biasadd_readvariableop_resource
identityИвconv1d/BiasAdd/ReadVariableOpв)conv1d/conv1d/ExpandDims_1/ReadVariableOpвconv1d_2/BiasAdd/ReadVariableOpв+conv1d_2/conv1d/ExpandDims_1/ReadVariableOpвdense/BiasAdd/ReadVariableOpвdense/MatMul/ReadVariableOpвdense_1/BiasAdd/ReadVariableOpвdense_1/MatMul/ReadVariableOpвdense_2/BiasAdd/ReadVariableOpвdense_2/MatMul/ReadVariableOpвdense_3/BiasAdd/ReadVariableOpвdense_3/MatMul/ReadVariableOpвdense_4/BiasAdd/ReadVariableOpвdense_4/MatMul/ReadVariableOpвdense_6/BiasAdd/ReadVariableOpвdense_6/MatMul/ReadVariableOpвdense_7/BiasAdd/ReadVariableOpвdense_7/MatMul/ReadVariableOpЗ
conv1d/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2
conv1d/conv1d/ExpandDims/dimм
conv1d/conv1d/ExpandDims
ExpandDimsinputs%conv1d/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d/conv1d/ExpandDims╧
)conv1d/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp2conv1d_conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02+
)conv1d/conv1d/ExpandDims_1/ReadVariableOpВ
conv1d/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2 
conv1d/conv1d/ExpandDims_1/dim╒
conv1d/conv1d/ExpandDims_1
ExpandDims1conv1d/conv1d/ExpandDims_1/ReadVariableOp:value:0'conv1d/conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d/conv1d/ExpandDims_1╘
conv1d/conv1dConv2D!conv1d/conv1d/ExpandDims:output:0#conv1d/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1d/conv1dи
conv1d/conv1d/SqueezeSqueezeconv1d/conv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d/conv1d/Squeezeв
conv1d/BiasAdd/ReadVariableOpReadVariableOp&conv1d_biasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02
conv1d/BiasAdd/ReadVariableOpй
conv1d/BiasAddBiasAddconv1d/conv1d/Squeeze:output:0%conv1d/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2
conv1d/BiasAddr
conv1d/ReluReluconv1d/BiasAdd:output:0*
T0*,
_output_shapes
:         а2
conv1d/Reluo
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
flatten/ConstФ
flatten/ReshapeReshapeconv1d/Relu:activations:0flatten/Const:output:0*
T0*)
_output_shapes
:         аГ2
flatten/Reshapeв
dense/MatMul/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
dense/MatMul/ReadVariableOpШ
dense/MatMulMatMulflatten/Reshape:output:0#dense/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense/MatMulЯ
dense/BiasAdd/ReadVariableOpReadVariableOp%dense_biasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02
dense/BiasAdd/ReadVariableOpЪ
dense/BiasAddBiasAdddense/MatMul:product:0$dense/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense/BiasAddk

dense/ReluReludense/BiasAdd:output:0*
T0*(
_output_shapes
:         Р2

dense/Reluж
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes
:	Р*
dtype02
dense_1/MatMul/ReadVariableOpЭ
dense_1/MatMulMatMuldense/Relu:activations:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
dense_1/MatMulд
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_1/BiasAdd/ReadVariableOpб
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
dense_1/BiasAddp
dense_1/TanhTanhdense_1/BiasAdd:output:0*
T0*'
_output_shapes
:         2
dense_1/Tanhе
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:<*
dtype02
dense_2/MatMul/ReadVariableOpХ
dense_2/MatMulMatMuldense_1/Tanh:y:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
dense_2/MatMulд
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:<*
dtype02 
dense_2/BiasAdd/ReadVariableOpб
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
dense_2/BiasAddp
dense_2/TanhTanhdense_2/BiasAdd:output:0*
T0*'
_output_shapes
:         <2
dense_2/Tanhж
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes
:	<░	*
dtype02
dense_3/MatMul/ReadVariableOpЦ
dense_3/MatMulMatMuldense_2/Tanh:y:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
dense_3/MatMulе
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes	
:░	*
dtype02 
dense_3/BiasAdd/ReadVariableOpв
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
dense_3/BiasAddq
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*(
_output_shapes
:         ░	2
dense_3/Reluз
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource* 
_output_shapes
:
░	╨A*
dtype02
dense_4/MatMul/ReadVariableOpа
dense_4/MatMulMatMuldense_3/Relu:activations:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
dense_4/MatMulе
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes	
:╨A*
dtype02 
dense_4/BiasAdd/ReadVariableOpв
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
dense_4/BiasAddq
dense_4/ReluReludense_4/BiasAdd:output:0*
T0*(
_output_shapes
:         ╨A2
dense_4/Relu
dropout/IdentityIdentitydense_4/Relu:activations:0*
T0*(
_output_shapes
:         ╨A2
dropout/Identityg
reshape/ShapeShapedropout/Identity:output:0*
T0*
_output_shapes
:2
reshape/ShapeД
reshape/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
reshape/strided_slice/stackИ
reshape/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
reshape/strided_slice/stack_1И
reshape/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
reshape/strided_slice/stack_2Т
reshape/strided_sliceStridedSlicereshape/Shape:output:0$reshape/strided_slice/stack:output:0&reshape/strided_slice/stack_1:output:0&reshape/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
reshape/strided_sliceu
reshape/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value
B :╨A2
reshape/Reshape/shape/1t
reshape/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
reshape/Reshape/shape/2╚
reshape/Reshape/shapePackreshape/strided_slice:output:0 reshape/Reshape/shape/1:output:0 reshape/Reshape/shape/2:output:0*
N*
T0*
_output_shapes
:2
reshape/Reshape/shapeЯ
reshape/ReshapeReshapedropout/Identity:output:0reshape/Reshape/shape:output:0*
T0*,
_output_shapes
:         ╨A2
reshape/ReshapeЛ
conv1d_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2 
conv1d_2/conv1d/ExpandDims/dim─
conv1d_2/conv1d/ExpandDims
ExpandDimsreshape/Reshape:output:0'conv1d_2/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2
conv1d_2/conv1d/ExpandDims╒
+conv1d_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp4conv1d_2_conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype02-
+conv1d_2/conv1d/ExpandDims_1/ReadVariableOpЖ
 conv1d_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2"
 conv1d_2/conv1d/ExpandDims_1/dim▌
conv1d_2/conv1d/ExpandDims_1
ExpandDims3conv1d_2/conv1d/ExpandDims_1/ReadVariableOp:value:0)conv1d_2/conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2
conv1d_2/conv1d/ExpandDims_1▄
conv1d_2/conv1dConv2D#conv1d_2/conv1d/ExpandDims:output:0%conv1d_2/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
conv1d_2/conv1dо
conv1d_2/conv1d/SqueezeSqueezeconv1d_2/conv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
conv1d_2/conv1d/Squeezeи
conv1d_2/BiasAdd/ReadVariableOpReadVariableOp(conv1d_2_biasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02!
conv1d_2/BiasAdd/ReadVariableOp▒
conv1d_2/BiasAddBiasAdd conv1d_2/conv1d/Squeeze:output:0'conv1d_2/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2
conv1d_2/BiasAddx
conv1d_2/ReluReluconv1d_2/BiasAdd:output:0*
T0*,
_output_shapes
:         а2
conv1d_2/Reluo
dense_5/ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
dense_5/ConstЦ
dense_5/ReshapeReshapeconv1d_2/Relu:activations:0dense_5/Const:output:0*
T0*)
_output_shapes
:         аГ2
dense_5/Reshapeи
dense_6/MatMul/ReadVariableOpReadVariableOp&dense_6_matmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02
dense_6/MatMul/ReadVariableOpЮ
dense_6/MatMulMatMuldense_5/Reshape:output:0%dense_6/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense_6/MatMulе
dense_6/BiasAdd/ReadVariableOpReadVariableOp'dense_6_biasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02 
dense_6/BiasAdd/ReadVariableOpв
dense_6/BiasAddBiasAdddense_6/MatMul:product:0&dense_6/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
dense_6/BiasAddq
dense_6/ReluReludense_6/BiasAdd:output:0*
T0*(
_output_shapes
:         Р2
dense_6/Reluж
dense_7/MatMul/ReadVariableOpReadVariableOp&dense_7_matmul_readvariableop_resource*
_output_shapes
:	РX*
dtype02
dense_7/MatMul/ReadVariableOpЯ
dense_7/MatMulMatMuldense_6/Relu:activations:0%dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
dense_7/MatMulд
dense_7/BiasAdd/ReadVariableOpReadVariableOp'dense_7_biasadd_readvariableop_resource*
_output_shapes
:X*
dtype02 
dense_7/BiasAdd/ReadVariableOpб
dense_7/BiasAddBiasAdddense_7/MatMul:product:0&dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
dense_7/BiasAddy
dense_7/SoftmaxSoftmaxdense_7/BiasAdd:output:0*
T0*'
_output_shapes
:         X2
dense_7/SoftmaxБ
dropout_2/IdentityIdentitydense_7/Softmax:softmax:0*
T0*'
_output_shapes
:         X2
dropout_2/Identity╬
IdentityIdentitydropout_2/Identity:output:0^conv1d/BiasAdd/ReadVariableOp*^conv1d/conv1d/ExpandDims_1/ReadVariableOp ^conv1d_2/BiasAdd/ReadVariableOp,^conv1d_2/conv1d/ExpandDims_1/ReadVariableOp^dense/BiasAdd/ReadVariableOp^dense/MatMul/ReadVariableOp^dense_1/BiasAdd/ReadVariableOp^dense_1/MatMul/ReadVariableOp^dense_2/BiasAdd/ReadVariableOp^dense_2/MatMul/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp^dense_6/BiasAdd/ReadVariableOp^dense_6/MatMul/ReadVariableOp^dense_7/BiasAdd/ReadVariableOp^dense_7/MatMul/ReadVariableOp*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::2>
conv1d/BiasAdd/ReadVariableOpconv1d/BiasAdd/ReadVariableOp2V
)conv1d/conv1d/ExpandDims_1/ReadVariableOp)conv1d/conv1d/ExpandDims_1/ReadVariableOp2B
conv1d_2/BiasAdd/ReadVariableOpconv1d_2/BiasAdd/ReadVariableOp2Z
+conv1d_2/conv1d/ExpandDims_1/ReadVariableOp+conv1d_2/conv1d/ExpandDims_1/ReadVariableOp2<
dense/BiasAdd/ReadVariableOpdense/BiasAdd/ReadVariableOp2:
dense/MatMul/ReadVariableOpdense/MatMul/ReadVariableOp2@
dense_1/BiasAdd/ReadVariableOpdense_1/BiasAdd/ReadVariableOp2>
dense_1/MatMul/ReadVariableOpdense_1/MatMul/ReadVariableOp2@
dense_2/BiasAdd/ReadVariableOpdense_2/BiasAdd/ReadVariableOp2>
dense_2/MatMul/ReadVariableOpdense_2/MatMul/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp2@
dense_6/BiasAdd/ReadVariableOpdense_6/BiasAdd/ReadVariableOp2>
dense_6/MatMul/ReadVariableOpdense_6/MatMul/ReadVariableOp2@
dense_7/BiasAdd/ReadVariableOpdense_7/BiasAdd/ReadVariableOp2>
dense_7/MatMul/ReadVariableOpdense_7/MatMul/ReadVariableOp:T P
,
_output_shapes
:         ╨A
 
_user_specified_nameinputs
▌
{
&__inference_dense_1_layer_call_fn_3235

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCallЎ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_1_layer_call_and_return_conditional_losses_23122
StatefulPartitionedCallО
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:         2

Identity"
identityIdentity:output:0*/
_input_shapes
:         Р::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         Р
 
_user_specified_nameinputs
Ж
`
A__inference_dropout_layer_call_and_return_conditional_losses_3307

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:         ╨A2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape╡
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:         ╨A*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *    2
dropout/GreaterEqual/y┐
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:         ╨A2
dropout/GreaterEqualА
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:         ╨A2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:         ╨A2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
╣?
Х
?__inference_model_layer_call_and_return_conditional_losses_2591
input_1
conv1d_2255
conv1d_2257

dense_2296

dense_2298
dense_1_2323
dense_1_2325
dense_2_2350
dense_2_2352
dense_3_2377
dense_3_2379
dense_4_2404
dense_4_2406
conv1d_2_2487
conv1d_2_2489
dense_6_2528
dense_6_2530
dense_7_2555
dense_7_2557
identityИвconv1d/StatefulPartitionedCallв conv1d_2/StatefulPartitionedCallвdense/StatefulPartitionedCallвdense_1/StatefulPartitionedCallвdense_2/StatefulPartitionedCallвdense_3/StatefulPartitionedCallвdense_4/StatefulPartitionedCallвdense_6/StatefulPartitionedCallвdense_7/StatefulPartitionedCallвdropout/StatefulPartitionedCallв!dropout_2/StatefulPartitionedCallП
conv1d/StatefulPartitionedCallStatefulPartitionedCallinput_1conv1d_2255conv1d_2257*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *I
fDRB
@__inference_conv1d_layer_call_and_return_conditional_losses_22442 
conv1d/StatefulPartitionedCallў
flatten/PartitionedCallPartitionedCall'conv1d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_flatten_layer_call_and_return_conditional_losses_22662
flatten/PartitionedCallЯ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
dense_2296
dense_2298*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_dense_layer_call_and_return_conditional_losses_22852
dense/StatefulPartitionedCallо
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_2323dense_1_2325*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_1_layer_call_and_return_conditional_losses_23122!
dense_1/StatefulPartitionedCall░
dense_2/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0dense_2_2350dense_2_2352*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         <*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_2_layer_call_and_return_conditional_losses_23392!
dense_2/StatefulPartitionedCall▒
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_2377dense_3_2379*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ░	*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_3_layer_call_and_return_conditional_losses_23662!
dense_3/StatefulPartitionedCall▒
dense_4/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0dense_4_2404dense_4_2406*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_4_layer_call_and_return_conditional_losses_23932!
dense_4/StatefulPartitionedCallП
dropout/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_24212!
dropout/StatefulPartitionedCall√
reshape/PartitionedCallPartitionedCall(dropout/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_reshape_layer_call_and_return_conditional_losses_24522
reshape/PartitionedCall▓
 conv1d_2/StatefulPartitionedCallStatefulPartitionedCall reshape/PartitionedCall:output:0conv1d_2_2487conv1d_2_2489*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         а*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *K
fFRD
B__inference_conv1d_2_layer_call_and_return_conditional_losses_24762"
 conv1d_2/StatefulPartitionedCall∙
dense_5/PartitionedCallPartitionedCall)conv1d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *)
_output_shapes
:         аГ* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_5_layer_call_and_return_conditional_losses_24982
dense_5/PartitionedCallй
dense_6/StatefulPartitionedCallStatefulPartitionedCall dense_5/PartitionedCall:output:0dense_6_2528dense_6_2530*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_6_layer_call_and_return_conditional_losses_25172!
dense_6/StatefulPartitionedCall░
dense_7/StatefulPartitionedCallStatefulPartitionedCall(dense_6/StatefulPartitionedCall:output:0dense_7_2555dense_7_2557*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_7_layer_call_and_return_conditional_losses_25442!
dense_7/StatefulPartitionedCall╢
!dropout_2/StatefulPartitionedCallStatefulPartitionedCall(dense_7/StatefulPartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         X* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *L
fGRE
C__inference_dropout_2_layer_call_and_return_conditional_losses_25722#
!dropout_2/StatefulPartitionedCallЇ
IdentityIdentity*dropout_2/StatefulPartitionedCall:output:0^conv1d/StatefulPartitionedCall!^conv1d_2/StatefulPartitionedCall^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall ^dense_2/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall ^dense_6/StatefulPartitionedCall ^dense_7/StatefulPartitionedCall ^dropout/StatefulPartitionedCall"^dropout_2/StatefulPartitionedCall*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::2@
conv1d/StatefulPartitionedCallconv1d/StatefulPartitionedCall2D
 conv1d_2/StatefulPartitionedCall conv1d_2/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2B
dense_6/StatefulPartitionedCalldense_6/StatefulPartitionedCall2B
dense_7/StatefulPartitionedCalldense_7/StatefulPartitionedCall2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_2/StatefulPartitionedCall!dropout_2/StatefulPartitionedCall:U Q
,
_output_shapes
:         ╨A
!
_user_specified_name	input_1
╡
]
A__inference_flatten_layer_call_and_return_conditional_losses_2266

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
Consti
ReshapeReshapeinputsConst:output:0*
T0*)
_output_shapes
:         аГ2	
Reshapef
IdentityIdentityReshape:output:0*
T0*)
_output_shapes
:         аГ2

Identity"
identityIdentity:output:0*+
_input_shapes
:         а:T P
,
_output_shapes
:         а
 
_user_specified_nameinputs
нВ
н
__inference__wrapped_model_2224
input_1<
8model_conv1d_conv1d_expanddims_1_readvariableop_resource0
,model_conv1d_biasadd_readvariableop_resource.
*model_dense_matmul_readvariableop_resource/
+model_dense_biasadd_readvariableop_resource0
,model_dense_1_matmul_readvariableop_resource1
-model_dense_1_biasadd_readvariableop_resource0
,model_dense_2_matmul_readvariableop_resource1
-model_dense_2_biasadd_readvariableop_resource0
,model_dense_3_matmul_readvariableop_resource1
-model_dense_3_biasadd_readvariableop_resource0
,model_dense_4_matmul_readvariableop_resource1
-model_dense_4_biasadd_readvariableop_resource>
:model_conv1d_2_conv1d_expanddims_1_readvariableop_resource2
.model_conv1d_2_biasadd_readvariableop_resource0
,model_dense_6_matmul_readvariableop_resource1
-model_dense_6_biasadd_readvariableop_resource0
,model_dense_7_matmul_readvariableop_resource1
-model_dense_7_biasadd_readvariableop_resource
identityИв#model/conv1d/BiasAdd/ReadVariableOpв/model/conv1d/conv1d/ExpandDims_1/ReadVariableOpв%model/conv1d_2/BiasAdd/ReadVariableOpв1model/conv1d_2/conv1d/ExpandDims_1/ReadVariableOpв"model/dense/BiasAdd/ReadVariableOpв!model/dense/MatMul/ReadVariableOpв$model/dense_1/BiasAdd/ReadVariableOpв#model/dense_1/MatMul/ReadVariableOpв$model/dense_2/BiasAdd/ReadVariableOpв#model/dense_2/MatMul/ReadVariableOpв$model/dense_3/BiasAdd/ReadVariableOpв#model/dense_3/MatMul/ReadVariableOpв$model/dense_4/BiasAdd/ReadVariableOpв#model/dense_4/MatMul/ReadVariableOpв$model/dense_6/BiasAdd/ReadVariableOpв#model/dense_6/MatMul/ReadVariableOpв$model/dense_7/BiasAdd/ReadVariableOpв#model/dense_7/MatMul/ReadVariableOpУ
"model/conv1d/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2$
"model/conv1d/conv1d/ExpandDims/dim┐
model/conv1d/conv1d/ExpandDims
ExpandDimsinput_1+model/conv1d/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2 
model/conv1d/conv1d/ExpandDimsс
/model/conv1d/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp8model_conv1d_conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype021
/model/conv1d/conv1d/ExpandDims_1/ReadVariableOpО
$model/conv1d/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2&
$model/conv1d/conv1d/ExpandDims_1/dimэ
 model/conv1d/conv1d/ExpandDims_1
ExpandDims7model/conv1d/conv1d/ExpandDims_1/ReadVariableOp:value:0-model/conv1d/conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2"
 model/conv1d/conv1d/ExpandDims_1ь
model/conv1d/conv1dConv2D'model/conv1d/conv1d/ExpandDims:output:0)model/conv1d/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
model/conv1d/conv1d║
model/conv1d/conv1d/SqueezeSqueezemodel/conv1d/conv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
model/conv1d/conv1d/Squeeze┤
#model/conv1d/BiasAdd/ReadVariableOpReadVariableOp,model_conv1d_biasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02%
#model/conv1d/BiasAdd/ReadVariableOp┴
model/conv1d/BiasAddBiasAdd$model/conv1d/conv1d/Squeeze:output:0+model/conv1d/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2
model/conv1d/BiasAddД
model/conv1d/ReluRelumodel/conv1d/BiasAdd:output:0*
T0*,
_output_shapes
:         а2
model/conv1d/Relu{
model/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
model/flatten/Constм
model/flatten/ReshapeReshapemodel/conv1d/Relu:activations:0model/flatten/Const:output:0*
T0*)
_output_shapes
:         аГ2
model/flatten/Reshape┤
!model/dense/MatMul/ReadVariableOpReadVariableOp*model_dense_matmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02#
!model/dense/MatMul/ReadVariableOp░
model/dense/MatMulMatMulmodel/flatten/Reshape:output:0)model/dense/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
model/dense/MatMul▒
"model/dense/BiasAdd/ReadVariableOpReadVariableOp+model_dense_biasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02$
"model/dense/BiasAdd/ReadVariableOp▓
model/dense/BiasAddBiasAddmodel/dense/MatMul:product:0*model/dense/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
model/dense/BiasAdd}
model/dense/ReluRelumodel/dense/BiasAdd:output:0*
T0*(
_output_shapes
:         Р2
model/dense/Relu╕
#model/dense_1/MatMul/ReadVariableOpReadVariableOp,model_dense_1_matmul_readvariableop_resource*
_output_shapes
:	Р*
dtype02%
#model/dense_1/MatMul/ReadVariableOp╡
model/dense_1/MatMulMatMulmodel/dense/Relu:activations:0+model/dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
model/dense_1/MatMul╢
$model/dense_1/BiasAdd/ReadVariableOpReadVariableOp-model_dense_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02&
$model/dense_1/BiasAdd/ReadVariableOp╣
model/dense_1/BiasAddBiasAddmodel/dense_1/MatMul:product:0,model/dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
model/dense_1/BiasAddВ
model/dense_1/TanhTanhmodel/dense_1/BiasAdd:output:0*
T0*'
_output_shapes
:         2
model/dense_1/Tanh╖
#model/dense_2/MatMul/ReadVariableOpReadVariableOp,model_dense_2_matmul_readvariableop_resource*
_output_shapes

:<*
dtype02%
#model/dense_2/MatMul/ReadVariableOpн
model/dense_2/MatMulMatMulmodel/dense_1/Tanh:y:0+model/dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
model/dense_2/MatMul╢
$model/dense_2/BiasAdd/ReadVariableOpReadVariableOp-model_dense_2_biasadd_readvariableop_resource*
_output_shapes
:<*
dtype02&
$model/dense_2/BiasAdd/ReadVariableOp╣
model/dense_2/BiasAddBiasAddmodel/dense_2/MatMul:product:0,model/dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         <2
model/dense_2/BiasAddВ
model/dense_2/TanhTanhmodel/dense_2/BiasAdd:output:0*
T0*'
_output_shapes
:         <2
model/dense_2/Tanh╕
#model/dense_3/MatMul/ReadVariableOpReadVariableOp,model_dense_3_matmul_readvariableop_resource*
_output_shapes
:	<░	*
dtype02%
#model/dense_3/MatMul/ReadVariableOpо
model/dense_3/MatMulMatMulmodel/dense_2/Tanh:y:0+model/dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
model/dense_3/MatMul╖
$model/dense_3/BiasAdd/ReadVariableOpReadVariableOp-model_dense_3_biasadd_readvariableop_resource*
_output_shapes	
:░	*
dtype02&
$model/dense_3/BiasAdd/ReadVariableOp║
model/dense_3/BiasAddBiasAddmodel/dense_3/MatMul:product:0,model/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
model/dense_3/BiasAddГ
model/dense_3/ReluRelumodel/dense_3/BiasAdd:output:0*
T0*(
_output_shapes
:         ░	2
model/dense_3/Relu╣
#model/dense_4/MatMul/ReadVariableOpReadVariableOp,model_dense_4_matmul_readvariableop_resource* 
_output_shapes
:
░	╨A*
dtype02%
#model/dense_4/MatMul/ReadVariableOp╕
model/dense_4/MatMulMatMul model/dense_3/Relu:activations:0+model/dense_4/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
model/dense_4/MatMul╖
$model/dense_4/BiasAdd/ReadVariableOpReadVariableOp-model_dense_4_biasadd_readvariableop_resource*
_output_shapes	
:╨A*
dtype02&
$model/dense_4/BiasAdd/ReadVariableOp║
model/dense_4/BiasAddBiasAddmodel/dense_4/MatMul:product:0,model/dense_4/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ╨A2
model/dense_4/BiasAddГ
model/dense_4/ReluRelumodel/dense_4/BiasAdd:output:0*
T0*(
_output_shapes
:         ╨A2
model/dense_4/ReluС
model/dropout/IdentityIdentity model/dense_4/Relu:activations:0*
T0*(
_output_shapes
:         ╨A2
model/dropout/Identityy
model/reshape/ShapeShapemodel/dropout/Identity:output:0*
T0*
_output_shapes
:2
model/reshape/ShapeР
!model/reshape/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2#
!model/reshape/strided_slice/stackФ
#model/reshape/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2%
#model/reshape/strided_slice/stack_1Ф
#model/reshape/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2%
#model/reshape/strided_slice/stack_2╢
model/reshape/strided_sliceStridedSlicemodel/reshape/Shape:output:0*model/reshape/strided_slice/stack:output:0,model/reshape/strided_slice/stack_1:output:0,model/reshape/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
model/reshape/strided_sliceБ
model/reshape/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value
B :╨A2
model/reshape/Reshape/shape/1А
model/reshape/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
model/reshape/Reshape/shape/2ц
model/reshape/Reshape/shapePack$model/reshape/strided_slice:output:0&model/reshape/Reshape/shape/1:output:0&model/reshape/Reshape/shape/2:output:0*
N*
T0*
_output_shapes
:2
model/reshape/Reshape/shape╖
model/reshape/ReshapeReshapemodel/dropout/Identity:output:0$model/reshape/Reshape/shape:output:0*
T0*,
_output_shapes
:         ╨A2
model/reshape/ReshapeЧ
$model/conv1d_2/conv1d/ExpandDims/dimConst*
_output_shapes
: *
dtype0*
valueB :
¤        2&
$model/conv1d_2/conv1d/ExpandDims/dim▄
 model/conv1d_2/conv1d/ExpandDims
ExpandDimsmodel/reshape/Reshape:output:0-model/conv1d_2/conv1d/ExpandDims/dim:output:0*
T0*0
_output_shapes
:         ╨A2"
 model/conv1d_2/conv1d/ExpandDimsч
1model/conv1d_2/conv1d/ExpandDims_1/ReadVariableOpReadVariableOp:model_conv1d_2_conv1d_expanddims_1_readvariableop_resource*$
_output_shapes
:░	а*
dtype023
1model/conv1d_2/conv1d/ExpandDims_1/ReadVariableOpТ
&model/conv1d_2/conv1d/ExpandDims_1/dimConst*
_output_shapes
: *
dtype0*
value	B : 2(
&model/conv1d_2/conv1d/ExpandDims_1/dimї
"model/conv1d_2/conv1d/ExpandDims_1
ExpandDims9model/conv1d_2/conv1d/ExpandDims_1/ReadVariableOp:value:0/model/conv1d_2/conv1d/ExpandDims_1/dim:output:0*
T0*(
_output_shapes
:░	а2$
"model/conv1d_2/conv1d/ExpandDims_1Ї
model/conv1d_2/conv1dConv2D)model/conv1d_2/conv1d/ExpandDims:output:0+model/conv1d_2/conv1d/ExpandDims_1:output:0*
T0*0
_output_shapes
:         а*
paddingSAME*
strides	
Р2
model/conv1d_2/conv1d└
model/conv1d_2/conv1d/SqueezeSqueezemodel/conv1d_2/conv1d:output:0*
T0*,
_output_shapes
:         а*
squeeze_dims

¤        2
model/conv1d_2/conv1d/Squeeze║
%model/conv1d_2/BiasAdd/ReadVariableOpReadVariableOp.model_conv1d_2_biasadd_readvariableop_resource*
_output_shapes	
:а*
dtype02'
%model/conv1d_2/BiasAdd/ReadVariableOp╔
model/conv1d_2/BiasAddBiasAdd&model/conv1d_2/conv1d/Squeeze:output:0-model/conv1d_2/BiasAdd/ReadVariableOp:value:0*
T0*,
_output_shapes
:         а2
model/conv1d_2/BiasAddК
model/conv1d_2/ReluRelumodel/conv1d_2/BiasAdd:output:0*
T0*,
_output_shapes
:         а2
model/conv1d_2/Relu{
model/dense_5/ConstConst*
_output_shapes
:*
dtype0*
valueB"    аA  2
model/dense_5/Constо
model/dense_5/ReshapeReshape!model/conv1d_2/Relu:activations:0model/dense_5/Const:output:0*
T0*)
_output_shapes
:         аГ2
model/dense_5/Reshape║
#model/dense_6/MatMul/ReadVariableOpReadVariableOp,model_dense_6_matmul_readvariableop_resource*!
_output_shapes
:аГР*
dtype02%
#model/dense_6/MatMul/ReadVariableOp╢
model/dense_6/MatMulMatMulmodel/dense_5/Reshape:output:0+model/dense_6/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
model/dense_6/MatMul╖
$model/dense_6/BiasAdd/ReadVariableOpReadVariableOp-model_dense_6_biasadd_readvariableop_resource*
_output_shapes	
:Р*
dtype02&
$model/dense_6/BiasAdd/ReadVariableOp║
model/dense_6/BiasAddBiasAddmodel/dense_6/MatMul:product:0,model/dense_6/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         Р2
model/dense_6/BiasAddГ
model/dense_6/ReluRelumodel/dense_6/BiasAdd:output:0*
T0*(
_output_shapes
:         Р2
model/dense_6/Relu╕
#model/dense_7/MatMul/ReadVariableOpReadVariableOp,model_dense_7_matmul_readvariableop_resource*
_output_shapes
:	РX*
dtype02%
#model/dense_7/MatMul/ReadVariableOp╖
model/dense_7/MatMulMatMul model/dense_6/Relu:activations:0+model/dense_7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
model/dense_7/MatMul╢
$model/dense_7/BiasAdd/ReadVariableOpReadVariableOp-model_dense_7_biasadd_readvariableop_resource*
_output_shapes
:X*
dtype02&
$model/dense_7/BiasAdd/ReadVariableOp╣
model/dense_7/BiasAddBiasAddmodel/dense_7/MatMul:product:0,model/dense_7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         X2
model/dense_7/BiasAddЛ
model/dense_7/SoftmaxSoftmaxmodel/dense_7/BiasAdd:output:0*
T0*'
_output_shapes
:         X2
model/dense_7/SoftmaxУ
model/dropout_2/IdentityIdentitymodel/dense_7/Softmax:softmax:0*
T0*'
_output_shapes
:         X2
model/dropout_2/Identity└
IdentityIdentity!model/dropout_2/Identity:output:0$^model/conv1d/BiasAdd/ReadVariableOp0^model/conv1d/conv1d/ExpandDims_1/ReadVariableOp&^model/conv1d_2/BiasAdd/ReadVariableOp2^model/conv1d_2/conv1d/ExpandDims_1/ReadVariableOp#^model/dense/BiasAdd/ReadVariableOp"^model/dense/MatMul/ReadVariableOp%^model/dense_1/BiasAdd/ReadVariableOp$^model/dense_1/MatMul/ReadVariableOp%^model/dense_2/BiasAdd/ReadVariableOp$^model/dense_2/MatMul/ReadVariableOp%^model/dense_3/BiasAdd/ReadVariableOp$^model/dense_3/MatMul/ReadVariableOp%^model/dense_4/BiasAdd/ReadVariableOp$^model/dense_4/MatMul/ReadVariableOp%^model/dense_6/BiasAdd/ReadVariableOp$^model/dense_6/MatMul/ReadVariableOp%^model/dense_7/BiasAdd/ReadVariableOp$^model/dense_7/MatMul/ReadVariableOp*
T0*'
_output_shapes
:         X2

Identity"
identityIdentity:output:0*s
_input_shapesb
`:         ╨A::::::::::::::::::2J
#model/conv1d/BiasAdd/ReadVariableOp#model/conv1d/BiasAdd/ReadVariableOp2b
/model/conv1d/conv1d/ExpandDims_1/ReadVariableOp/model/conv1d/conv1d/ExpandDims_1/ReadVariableOp2N
%model/conv1d_2/BiasAdd/ReadVariableOp%model/conv1d_2/BiasAdd/ReadVariableOp2f
1model/conv1d_2/conv1d/ExpandDims_1/ReadVariableOp1model/conv1d_2/conv1d/ExpandDims_1/ReadVariableOp2H
"model/dense/BiasAdd/ReadVariableOp"model/dense/BiasAdd/ReadVariableOp2F
!model/dense/MatMul/ReadVariableOp!model/dense/MatMul/ReadVariableOp2L
$model/dense_1/BiasAdd/ReadVariableOp$model/dense_1/BiasAdd/ReadVariableOp2J
#model/dense_1/MatMul/ReadVariableOp#model/dense_1/MatMul/ReadVariableOp2L
$model/dense_2/BiasAdd/ReadVariableOp$model/dense_2/BiasAdd/ReadVariableOp2J
#model/dense_2/MatMul/ReadVariableOp#model/dense_2/MatMul/ReadVariableOp2L
$model/dense_3/BiasAdd/ReadVariableOp$model/dense_3/BiasAdd/ReadVariableOp2J
#model/dense_3/MatMul/ReadVariableOp#model/dense_3/MatMul/ReadVariableOp2L
$model/dense_4/BiasAdd/ReadVariableOp$model/dense_4/BiasAdd/ReadVariableOp2J
#model/dense_4/MatMul/ReadVariableOp#model/dense_4/MatMul/ReadVariableOp2L
$model/dense_6/BiasAdd/ReadVariableOp$model/dense_6/BiasAdd/ReadVariableOp2J
#model/dense_6/MatMul/ReadVariableOp#model/dense_6/MatMul/ReadVariableOp2L
$model/dense_7/BiasAdd/ReadVariableOp$model/dense_7/BiasAdd/ReadVariableOp2J
#model/dense_7/MatMul/ReadVariableOp#model/dense_7/MatMul/ReadVariableOp:U Q
,
_output_shapes
:         ╨A
!
_user_specified_name	input_1
б
_
&__inference_dropout_layer_call_fn_3317

inputs
identityИвStatefulPartitionedCall▌
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_24212
StatefulPartitionedCallП
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*'
_input_shapes
:         ╨A22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
▀
{
&__inference_dense_4_layer_call_fn_3295

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCallў
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ╨A*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_dense_4_layer_call_and_return_conditional_losses_23932
StatefulPartitionedCallП
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*/
_input_shapes
:         ░	::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         ░	
 
_user_specified_nameinputs
ё	
┌
A__inference_dense_3_layer_call_and_return_conditional_losses_3266

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	<░	*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:░	*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ░	2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         ░	2
ReluШ
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:         ░	2

Identity"
identityIdentity:output:0*.
_input_shapes
:         <::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         <
 
_user_specified_nameinputs
Э
B
&__inference_reshape_layer_call_fn_3340

inputs
identity╔
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:         ╨A* 
_read_only_resource_inputs
 *2
config_proto" 

CPU

GPU2*0,1J 8В *J
fERC
A__inference_reshape_layer_call_and_return_conditional_losses_24522
PartitionedCallq
IdentityIdentityPartitionedCall:output:0*
T0*,
_output_shapes
:         ╨A2

Identity"
identityIdentity:output:0*'
_input_shapes
:         ╨A:P L
(
_output_shapes
:         ╨A
 
_user_specified_nameinputs
╣Н
┼
 __inference__traced_restore_3680
file_prefix"
assignvariableop_conv1d_kernel"
assignvariableop_1_conv1d_bias#
assignvariableop_2_dense_kernel!
assignvariableop_3_dense_bias%
!assignvariableop_4_dense_1_kernel#
assignvariableop_5_dense_1_bias%
!assignvariableop_6_dense_2_kernel#
assignvariableop_7_dense_2_bias%
!assignvariableop_8_dense_3_kernel#
assignvariableop_9_dense_3_bias&
"assignvariableop_10_dense_4_kernel$
 assignvariableop_11_dense_4_bias'
#assignvariableop_12_conv1d_2_kernel%
!assignvariableop_13_conv1d_2_bias&
"assignvariableop_14_dense_6_kernel$
 assignvariableop_15_dense_6_bias&
"assignvariableop_16_dense_7_kernel$
 assignvariableop_17_dense_7_bias
assignvariableop_18_total
assignvariableop_19_count
assignvariableop_20_total_1
assignvariableop_21_count_1)
%assignvariableop_22_conv1d_2_kernel_m'
#assignvariableop_23_conv1d_2_bias_m(
$assignvariableop_24_dense_6_kernel_m&
"assignvariableop_25_dense_6_bias_m(
$assignvariableop_26_dense_7_kernel_m&
"assignvariableop_27_dense_7_bias_m)
%assignvariableop_28_conv1d_2_kernel_v'
#assignvariableop_29_conv1d_2_bias_v(
$assignvariableop_30_dense_6_kernel_v&
"assignvariableop_31_dense_6_bias_v(
$assignvariableop_32_dense_7_kernel_v&
"assignvariableop_33_dense_7_bias_v
identity_35ИвAssignVariableOpвAssignVariableOp_1вAssignVariableOp_10вAssignVariableOp_11вAssignVariableOp_12вAssignVariableOp_13вAssignVariableOp_14вAssignVariableOp_15вAssignVariableOp_16вAssignVariableOp_17вAssignVariableOp_18вAssignVariableOp_19вAssignVariableOp_2вAssignVariableOp_20вAssignVariableOp_21вAssignVariableOp_22вAssignVariableOp_23вAssignVariableOp_24вAssignVariableOp_25вAssignVariableOp_26вAssignVariableOp_27вAssignVariableOp_28вAssignVariableOp_29вAssignVariableOp_3вAssignVariableOp_30вAssignVariableOp_31вAssignVariableOp_32вAssignVariableOp_33вAssignVariableOp_4вAssignVariableOp_5вAssignVariableOp_6вAssignVariableOp_7вAssignVariableOp_8вAssignVariableOp_9┴
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:#*
dtype0*═
value├B└#B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_names╘
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:#*
dtype0*Y
valuePBN#B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices▌
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*в
_output_shapesП
М:::::::::::::::::::::::::::::::::::*1
dtypes'
%2#2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

IdentityЭ
AssignVariableOpAssignVariableOpassignvariableop_conv1d_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1г
AssignVariableOp_1AssignVariableOpassignvariableop_1_conv1d_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2д
AssignVariableOp_2AssignVariableOpassignvariableop_2_dense_kernelIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3в
AssignVariableOp_3AssignVariableOpassignvariableop_3_dense_biasIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4ж
AssignVariableOp_4AssignVariableOp!assignvariableop_4_dense_1_kernelIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5д
AssignVariableOp_5AssignVariableOpassignvariableop_5_dense_1_biasIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6ж
AssignVariableOp_6AssignVariableOp!assignvariableop_6_dense_2_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7д
AssignVariableOp_7AssignVariableOpassignvariableop_7_dense_2_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8ж
AssignVariableOp_8AssignVariableOp!assignvariableop_8_dense_3_kernelIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9д
AssignVariableOp_9AssignVariableOpassignvariableop_9_dense_3_biasIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10к
AssignVariableOp_10AssignVariableOp"assignvariableop_10_dense_4_kernelIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11и
AssignVariableOp_11AssignVariableOp assignvariableop_11_dense_4_biasIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:2
Identity_12л
AssignVariableOp_12AssignVariableOp#assignvariableop_12_conv1d_2_kernelIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13й
AssignVariableOp_13AssignVariableOp!assignvariableop_13_conv1d_2_biasIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14к
AssignVariableOp_14AssignVariableOp"assignvariableop_14_dense_6_kernelIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15и
AssignVariableOp_15AssignVariableOp assignvariableop_15_dense_6_biasIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16к
AssignVariableOp_16AssignVariableOp"assignvariableop_16_dense_7_kernelIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17и
AssignVariableOp_17AssignVariableOp assignvariableop_17_dense_7_biasIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18б
AssignVariableOp_18AssignVariableOpassignvariableop_18_totalIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19б
AssignVariableOp_19AssignVariableOpassignvariableop_19_countIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:2
Identity_20г
AssignVariableOp_20AssignVariableOpassignvariableop_20_total_1Identity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21г
AssignVariableOp_21AssignVariableOpassignvariableop_21_count_1Identity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22н
AssignVariableOp_22AssignVariableOp%assignvariableop_22_conv1d_2_kernel_mIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23л
AssignVariableOp_23AssignVariableOp#assignvariableop_23_conv1d_2_bias_mIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24м
AssignVariableOp_24AssignVariableOp$assignvariableop_24_dense_6_kernel_mIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25к
AssignVariableOp_25AssignVariableOp"assignvariableop_25_dense_6_bias_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26м
AssignVariableOp_26AssignVariableOp$assignvariableop_26_dense_7_kernel_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27к
AssignVariableOp_27AssignVariableOp"assignvariableop_27_dense_7_bias_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28н
AssignVariableOp_28AssignVariableOp%assignvariableop_28_conv1d_2_kernel_vIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29л
AssignVariableOp_29AssignVariableOp#assignvariableop_29_conv1d_2_bias_vIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30м
AssignVariableOp_30AssignVariableOp$assignvariableop_30_dense_6_kernel_vIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31к
AssignVariableOp_31AssignVariableOp"assignvariableop_31_dense_6_bias_vIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32м
AssignVariableOp_32AssignVariableOp$assignvariableop_32_dense_7_kernel_vIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33к
AssignVariableOp_33AssignVariableOp"assignvariableop_33_dense_7_bias_vIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_339
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp╩
Identity_34Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_34╜
Identity_35IdentityIdentity_34:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
T0*
_output_shapes
: 2
Identity_35"#
identity_35Identity_35:output:0*Я
_input_shapesН
К: ::::::::::::::::::::::::::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
▌
y
$__inference_dense_layer_call_fn_3215

inputs
unknown
	unknown_0
identityИвStatefulPartitionedCallї
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         Р*$
_read_only_resource_inputs
*2
config_proto" 

CPU

GPU2*0,1J 8В *H
fCRA
?__inference_dense_layer_call_and_return_conditional_losses_22852
StatefulPartitionedCallП
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:         Р2

Identity"
identityIdentity:output:0*0
_input_shapes
:         аГ::22
StatefulPartitionedCallStatefulPartitionedCall:Q M
)
_output_shapes
:         аГ
 
_user_specified_nameinputs"▒L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*▒
serving_defaultЭ
@
input_15
serving_default_input_1:0         ╨A=
	dropout_20
StatefulPartitionedCall:0         Xtensorflow/serving/predict:зи
╕r
layer-0
layer_with_weights-0
layer-1
layer-2
layer_with_weights-1
layer-3
layer_with_weights-2
layer-4
layer_with_weights-3
layer-5
layer_with_weights-4
layer-6
layer_with_weights-5
layer-7
	layer-8

layer-9
layer_with_weights-6
layer-10
layer-11
layer_with_weights-7
layer-12
layer_with_weights-8
layer-13
layer-14
	optimizer

signatures
#_self_saveable_object_factories
trainable_variables
regularization_losses
	variables
	keras_api
╥_default_save_signature
+╙&call_and_return_all_conditional_losses
╘__call__"Уm
_tf_keras_networkўl{"class_name": "Functional", "name": "model", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "model", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv1D", "config": {"name": "conv1d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "filters": 800, "kernel_size": {"class_name": "__tuple__", "items": [1200]}, "strides": {"class_name": "__tuple__", "items": [400]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv1d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["conv1d", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 20, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_1", "inbound_nodes": [[["dense", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 60, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_2", "inbound_nodes": [[["dense_1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 1200, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_3", "inbound_nodes": [[["dense_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 8400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_4", "inbound_nodes": [[["dense_3", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.0, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["dense_4", 0, 0, {}]]]}, {"class_name": "Reshape", "config": {"name": "reshape", "trainable": true, "dtype": "float32", "target_shape": {"class_name": "__tuple__", "items": [8400, 1]}}, "name": "reshape", "inbound_nodes": [[["dropout", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "conv1d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "filters": 800, "kernel_size": {"class_name": "__tuple__", "items": [1200]}, "strides": {"class_name": "__tuple__", "items": [400]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv1d_2", "inbound_nodes": [[["reshape", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "dense_5", "inbound_nodes": [[["conv1d_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_6", "trainable": true, "dtype": "float32", "units": 400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_6", "inbound_nodes": [[["dense_5", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_7", "trainable": true, "dtype": "float32", "units": 88, "activation": "softmax", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_7", "inbound_nodes": [[["dense_6", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.0, "noise_shape": null, "seed": null}, "name": "dropout_2", "inbound_nodes": [[["dense_7", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["dropout_2", 0, 0]]}, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 8400, 1]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "model", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv1D", "config": {"name": "conv1d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "filters": 800, "kernel_size": {"class_name": "__tuple__", "items": [1200]}, "strides": {"class_name": "__tuple__", "items": [400]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv1d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["conv1d", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 20, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_1", "inbound_nodes": [[["dense", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 60, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_2", "inbound_nodes": [[["dense_1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 1200, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_3", "inbound_nodes": [[["dense_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 8400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_4", "inbound_nodes": [[["dense_3", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.0, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["dense_4", 0, 0, {}]]]}, {"class_name": "Reshape", "config": {"name": "reshape", "trainable": true, "dtype": "float32", "target_shape": {"class_name": "__tuple__", "items": [8400, 1]}}, "name": "reshape", "inbound_nodes": [[["dropout", 0, 0, {}]]]}, {"class_name": "Conv1D", "config": {"name": "conv1d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "filters": 800, "kernel_size": {"class_name": "__tuple__", "items": [1200]}, "strides": {"class_name": "__tuple__", "items": [400]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv1d_2", "inbound_nodes": [[["reshape", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "dense_5", "inbound_nodes": [[["conv1d_2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_6", "trainable": true, "dtype": "float32", "units": 400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_6", "inbound_nodes": [[["dense_5", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_7", "trainable": true, "dtype": "float32", "units": 88, "activation": "softmax", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_7", "inbound_nodes": [[["dense_6", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.0, "noise_shape": null, "seed": null}, "name": "dropout_2", "inbound_nodes": [[["dense_7", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["dropout_2", 0, 0]]}}, "training_config": {"loss": {"class_name": "CategoricalCrossentropy", "config": {"reduction": "auto", "name": "categorical_crossentropy", "from_logits": false, "label_smoothing": 0}}, "metrics": [[{"class_name": "CategoricalAccuracy", "config": {"name": "categorical_accuracy", "dtype": "float32"}}]], "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adamax", "config": {"name": "Adamax", "learning_rate": 9.999999747378752e-05, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07}}}}
Ъ
#_self_saveable_object_factories"Є
_tf_keras_input_layer╥{"class_name": "InputLayer", "name": "input_1", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}}
Л

kernel
bias
#_self_saveable_object_factories
trainable_variables
regularization_losses
	variables
	keras_api
+╒&call_and_return_all_conditional_losses
╓__call__"┐	
_tf_keras_layerе	{"class_name": "Conv1D", "name": "conv1d", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv1d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "filters": 800, "kernel_size": {"class_name": "__tuple__", "items": [1200]}, "strides": {"class_name": "__tuple__", "items": [400]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 1}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 8400, 1]}}
Й
#_self_saveable_object_factories
 trainable_variables
!regularization_losses
"	variables
#	keras_api
+╫&call_and_return_all_conditional_losses
╪__call__"╙
_tf_keras_layer╣{"class_name": "Flatten", "name": "flatten", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}}}
Ъ

$kernel
%bias
#&_self_saveable_object_factories
'trainable_variables
(regularization_losses
)	variables
*	keras_api
+┘&call_and_return_all_conditional_losses
┌__call__"╬
_tf_keras_layer┤{"class_name": "Dense", "name": "dense", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 16800}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 16800]}}
Щ

+kernel
,bias
#-_self_saveable_object_factories
.trainable_variables
/regularization_losses
0	variables
1	keras_api
+█&call_and_return_all_conditional_losses
▄__call__"═
_tf_keras_layer│{"class_name": "Dense", "name": "dense_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 20, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 400}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 400]}}
Ч

2kernel
3bias
#4_self_saveable_object_factories
5trainable_variables
6regularization_losses
7	variables
8	keras_api
+▌&call_and_return_all_conditional_losses
▐__call__"╦
_tf_keras_layer▒{"class_name": "Dense", "name": "dense_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 60, "activation": "tanh", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 20]}}
Щ

9kernel
:bias
#;_self_saveable_object_factories
<trainable_variables
=regularization_losses
>	variables
?	keras_api
+▀&call_and_return_all_conditional_losses
р__call__"═
_tf_keras_layer│{"class_name": "Dense", "name": "dense_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 1200, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 60}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 60]}}
Э

@kernel
Abias
#B_self_saveable_object_factories
Ctrainable_variables
Dregularization_losses
E	variables
F	keras_api
+с&call_and_return_all_conditional_losses
т__call__"╤
_tf_keras_layer╖{"class_name": "Dense", "name": "dense_4", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 8400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1200}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 1200]}}
И
#G_self_saveable_object_factories
Htrainable_variables
Iregularization_losses
J	variables
K	keras_api
+у&call_and_return_all_conditional_losses
ф__call__"╥
_tf_keras_layer╕{"class_name": "Dropout", "name": "dropout", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.0, "noise_shape": null, "seed": null}}
Ъ
#L_self_saveable_object_factories
Mtrainable_variables
Nregularization_losses
O	variables
P	keras_api
+х&call_and_return_all_conditional_losses
ц__call__"ф
_tf_keras_layer╩{"class_name": "Reshape", "name": "reshape", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "reshape", "trainable": true, "dtype": "float32", "target_shape": {"class_name": "__tuple__", "items": [8400, 1]}}}
П

Qkernel
Rbias
#S_self_saveable_object_factories
Ttrainable_variables
Uregularization_losses
V	variables
W	keras_api
+ч&call_and_return_all_conditional_losses
ш__call__"├	
_tf_keras_layerй	{"class_name": "Conv1D", "name": "conv1d_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv1d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 8400, 1]}, "dtype": "float32", "filters": 800, "kernel_size": {"class_name": "__tuple__", "items": [1200]}, "strides": {"class_name": "__tuple__", "items": [400]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 3, "axes": {"-1": 1}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 8400, 1]}}
Й
#X_self_saveable_object_factories
Ytrainable_variables
Zregularization_losses
[	variables
\	keras_api
+щ&call_and_return_all_conditional_losses
ъ__call__"╙
_tf_keras_layer╣{"class_name": "Flatten", "name": "dense_5", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}}}
Ю

]kernel
^bias
#__self_saveable_object_factories
`trainable_variables
aregularization_losses
b	variables
c	keras_api
+ы&call_and_return_all_conditional_losses
ь__call__"╥
_tf_keras_layer╕{"class_name": "Dense", "name": "dense_6", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_6", "trainable": true, "dtype": "float32", "units": 400, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 16800}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 16800]}}
Ь

dkernel
ebias
#f_self_saveable_object_factories
gtrainable_variables
hregularization_losses
i	variables
j	keras_api
+э&call_and_return_all_conditional_losses
ю__call__"╨
_tf_keras_layer╢{"class_name": "Dense", "name": "dense_7", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_7", "trainable": true, "dtype": "float32", "units": 88, "activation": "softmax", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 400}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 400]}}
М
#k_self_saveable_object_factories
ltrainable_variables
mregularization_losses
n	variables
o	keras_api
+я&call_and_return_all_conditional_losses
Ё__call__"╓
_tf_keras_layer╝{"class_name": "Dropout", "name": "dropout_2", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dropout_2", "trainable": true, "dtype": "float32", "rate": 0.0, "noise_shape": null, "seed": null}}
ЛQm╞Rm╟]m╚^m╔dm╩em╦Qv╠Rv═]v╬^v╧dv╨ev╤"
	optimizer
-
ёserving_default"
signature_map
 "
trackable_dict_wrapper
ж
0
1
$2
%3
+4
,5
26
37
98
:9
@10
A11
Q12
R13
]14
^15
d16
e17"
trackable_list_wrapper
 "
trackable_list_wrapper
ж
0
1
$2
%3
+4
,5
26
37
98
:9
@10
A11
Q12
R13
]14
^15
d16
e17"
trackable_list_wrapper
╬
player_metrics
trainable_variables
qmetrics
regularization_losses
rnon_trainable_variables

slayers
	variables
tlayer_regularization_losses
╘__call__
╥_default_save_signature
+╙&call_and_return_all_conditional_losses
'╙"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
%:#░	а2conv1d/kernel
:а2conv1d/bias
 "
trackable_dict_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
░
ulayer_metrics
trainable_variables
vmetrics
regularization_losses
wnon_trainable_variables

xlayers
	variables
ylayer_regularization_losses
╓__call__
+╒&call_and_return_all_conditional_losses
'╒"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
░
zlayer_metrics
 trainable_variables
{metrics
!regularization_losses
|non_trainable_variables

}layers
"	variables
~layer_regularization_losses
╪__call__
+╫&call_and_return_all_conditional_losses
'╫"call_and_return_conditional_losses"
_generic_user_object
!:аГР2dense/kernel
:Р2
dense/bias
 "
trackable_dict_wrapper
.
$0
%1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
$0
%1"
trackable_list_wrapper
┤
layer_metrics
'trainable_variables
Аmetrics
(regularization_losses
Бnon_trainable_variables
Вlayers
)	variables
 Гlayer_regularization_losses
┌__call__
+┘&call_and_return_all_conditional_losses
'┘"call_and_return_conditional_losses"
_generic_user_object
!:	Р2dense_1/kernel
:2dense_1/bias
 "
trackable_dict_wrapper
.
+0
,1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
╡
Дlayer_metrics
.trainable_variables
Еmetrics
/regularization_losses
Жnon_trainable_variables
Зlayers
0	variables
 Иlayer_regularization_losses
▄__call__
+█&call_and_return_all_conditional_losses
'█"call_and_return_conditional_losses"
_generic_user_object
 :<2dense_2/kernel
:<2dense_2/bias
 "
trackable_dict_wrapper
.
20
31"
trackable_list_wrapper
 "
trackable_list_wrapper
.
20
31"
trackable_list_wrapper
╡
Йlayer_metrics
5trainable_variables
Кmetrics
6regularization_losses
Лnon_trainable_variables
Мlayers
7	variables
 Нlayer_regularization_losses
▐__call__
+▌&call_and_return_all_conditional_losses
'▌"call_and_return_conditional_losses"
_generic_user_object
!:	<░	2dense_3/kernel
:░	2dense_3/bias
 "
trackable_dict_wrapper
.
90
:1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
90
:1"
trackable_list_wrapper
╡
Оlayer_metrics
<trainable_variables
Пmetrics
=regularization_losses
Рnon_trainable_variables
Сlayers
>	variables
 Тlayer_regularization_losses
р__call__
+▀&call_and_return_all_conditional_losses
'▀"call_and_return_conditional_losses"
_generic_user_object
": 
░	╨A2dense_4/kernel
:╨A2dense_4/bias
 "
trackable_dict_wrapper
.
@0
A1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
@0
A1"
trackable_list_wrapper
╡
Уlayer_metrics
Ctrainable_variables
Фmetrics
Dregularization_losses
Хnon_trainable_variables
Цlayers
E	variables
 Чlayer_regularization_losses
т__call__
+с&call_and_return_all_conditional_losses
'с"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
Шlayer_metrics
Htrainable_variables
Щmetrics
Iregularization_losses
Ъnon_trainable_variables
Ыlayers
J	variables
 Ьlayer_regularization_losses
ф__call__
+у&call_and_return_all_conditional_losses
'у"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
Эlayer_metrics
Mtrainable_variables
Юmetrics
Nregularization_losses
Яnon_trainable_variables
аlayers
O	variables
 бlayer_regularization_losses
ц__call__
+х&call_and_return_all_conditional_losses
'х"call_and_return_conditional_losses"
_generic_user_object
':%░	а2conv1d_2/kernel
:а2conv1d_2/bias
 "
trackable_dict_wrapper
.
Q0
R1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
Q0
R1"
trackable_list_wrapper
╡
вlayer_metrics
Ttrainable_variables
гmetrics
Uregularization_losses
дnon_trainable_variables
еlayers
V	variables
 жlayer_regularization_losses
ш__call__
+ч&call_and_return_all_conditional_losses
'ч"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
зlayer_metrics
Ytrainable_variables
иmetrics
Zregularization_losses
йnon_trainable_variables
кlayers
[	variables
 лlayer_regularization_losses
ъ__call__
+щ&call_and_return_all_conditional_losses
'щ"call_and_return_conditional_losses"
_generic_user_object
#:!аГР2dense_6/kernel
:Р2dense_6/bias
 "
trackable_dict_wrapper
.
]0
^1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
]0
^1"
trackable_list_wrapper
╡
мlayer_metrics
`trainable_variables
нmetrics
aregularization_losses
оnon_trainable_variables
пlayers
b	variables
 ░layer_regularization_losses
ь__call__
+ы&call_and_return_all_conditional_losses
'ы"call_and_return_conditional_losses"
_generic_user_object
!:	РX2dense_7/kernel
:X2dense_7/bias
 "
trackable_dict_wrapper
.
d0
e1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
d0
e1"
trackable_list_wrapper
╡
▒layer_metrics
gtrainable_variables
▓metrics
hregularization_losses
│non_trainable_variables
┤layers
i	variables
 ╡layer_regularization_losses
ю__call__
+э&call_and_return_all_conditional_losses
'э"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
╢layer_metrics
ltrainable_variables
╖metrics
mregularization_losses
╕non_trainable_variables
╣layers
n	variables
 ║layer_regularization_losses
Ё__call__
+я&call_and_return_all_conditional_losses
'я"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
0
╗0
╝1"
trackable_list_wrapper
 "
trackable_list_wrapper
О
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
┐

╜total

╛count
┐	variables
└	keras_api"Д
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
А

┴total

┬count
├
_fn_kwargs
─	variables
┼	keras_api"┤
_tf_keras_metricЩ{"class_name": "CategoricalAccuracy", "name": "categorical_accuracy", "dtype": "float32", "config": {"name": "categorical_accuracy", "dtype": "float32"}}
:  (2total
:  (2count
0
╜0
╛1"
trackable_list_wrapper
.
┐	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
┴0
┬1"
trackable_list_wrapper
.
─	variables"
_generic_user_object
':%░	а2conv1d_2/kernel/m
:а2conv1d_2/bias/m
#:!аГР2dense_6/kernel/m
:Р2dense_6/bias/m
!:	РX2dense_7/kernel/m
:X2dense_7/bias/m
':%░	а2conv1d_2/kernel/v
:а2conv1d_2/bias/v
#:!аГР2dense_6/kernel/v
:Р2dense_6/bias/v
!:	РX2dense_7/kernel/v
:X2dense_7/bias/v
т2▀
__inference__wrapped_model_2224╗
Л▓З
FullArgSpec
argsЪ 
varargsjargs
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *+в(
&К#
input_1         ╨A
╩2╟
?__inference_model_layer_call_and_return_conditional_losses_3077
?__inference_model_layer_call_and_return_conditional_losses_2985
?__inference_model_layer_call_and_return_conditional_losses_2591
?__inference_model_layer_call_and_return_conditional_losses_2645└
╖▓│
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
▐2█
$__inference_model_layer_call_fn_2836
$__inference_model_layer_call_fn_2741
$__inference_model_layer_call_fn_3159
$__inference_model_layer_call_fn_3118└
╖▓│
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
ъ2ч
@__inference_conv1d_layer_call_and_return_conditional_losses_3175в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╧2╠
%__inference_conv1d_layer_call_fn_3184в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_flatten_layer_call_and_return_conditional_losses_3190в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_flatten_layer_call_fn_3195в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
щ2ц
?__inference_dense_layer_call_and_return_conditional_losses_3206в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╬2╦
$__inference_dense_layer_call_fn_3215в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_dense_1_layer_call_and_return_conditional_losses_3226в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_dense_1_layer_call_fn_3235в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_dense_2_layer_call_and_return_conditional_losses_3246в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_dense_2_layer_call_fn_3255в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_dense_3_layer_call_and_return_conditional_losses_3266в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_dense_3_layer_call_fn_3275в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_dense_4_layer_call_and_return_conditional_losses_3286в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_dense_4_layer_call_fn_3295в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
└2╜
A__inference_dropout_layer_call_and_return_conditional_losses_3307
A__inference_dropout_layer_call_and_return_conditional_losses_3312┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
К2З
&__inference_dropout_layer_call_fn_3317
&__inference_dropout_layer_call_fn_3322┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
ы2ш
A__inference_reshape_layer_call_and_return_conditional_losses_3335в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_reshape_layer_call_fn_3340в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ь2щ
B__inference_conv1d_2_layer_call_and_return_conditional_losses_3356в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╤2╬
'__inference_conv1d_2_layer_call_fn_3365в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_dense_5_layer_call_and_return_conditional_losses_3371в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_dense_5_layer_call_fn_3376в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_dense_6_layer_call_and_return_conditional_losses_3387в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_dense_6_layer_call_fn_3396в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_dense_7_layer_call_and_return_conditional_losses_3407в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_dense_7_layer_call_fn_3416в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
─2┴
C__inference_dropout_2_layer_call_and_return_conditional_losses_3428
C__inference_dropout_2_layer_call_and_return_conditional_losses_3433┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
О2Л
(__inference_dropout_2_layer_call_fn_3443
(__inference_dropout_2_layer_call_fn_3438┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╔B╞
"__inference_signature_wrapper_2879input_1"Ф
Н▓Й
FullArgSpec
argsЪ 
varargs
 
varkwjkwargs
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 ж
__inference__wrapped_model_2224В$%+,239:@AQR]^de5в2
+в(
&К#
input_1         ╨A
к "5к2
0
	dropout_2#К 
	dropout_2         Xм
B__inference_conv1d_2_layer_call_and_return_conditional_losses_3356fQR4в1
*в'
%К"
inputs         ╨A
к "*в'
 К
0         а
Ъ Д
'__inference_conv1d_2_layer_call_fn_3365YQR4в1
*в'
%К"
inputs         ╨A
к "К         ак
@__inference_conv1d_layer_call_and_return_conditional_losses_3175f4в1
*в'
%К"
inputs         ╨A
к "*в'
 К
0         а
Ъ В
%__inference_conv1d_layer_call_fn_3184Y4в1
*в'
%К"
inputs         ╨A
к "К         ав
A__inference_dense_1_layer_call_and_return_conditional_losses_3226]+,0в-
&в#
!К
inputs         Р
к "%в"
К
0         
Ъ z
&__inference_dense_1_layer_call_fn_3235P+,0в-
&в#
!К
inputs         Р
к "К         б
A__inference_dense_2_layer_call_and_return_conditional_losses_3246\23/в,
%в"
 К
inputs         
к "%в"
К
0         <
Ъ y
&__inference_dense_2_layer_call_fn_3255O23/в,
%в"
 К
inputs         
к "К         <в
A__inference_dense_3_layer_call_and_return_conditional_losses_3266]9:/в,
%в"
 К
inputs         <
к "&в#
К
0         ░	
Ъ z
&__inference_dense_3_layer_call_fn_3275P9:/в,
%в"
 К
inputs         <
к "К         ░	г
A__inference_dense_4_layer_call_and_return_conditional_losses_3286^@A0в-
&в#
!К
inputs         ░	
к "&в#
К
0         ╨A
Ъ {
&__inference_dense_4_layer_call_fn_3295Q@A0в-
&в#
!К
inputs         ░	
к "К         ╨Aд
A__inference_dense_5_layer_call_and_return_conditional_losses_3371_4в1
*в'
%К"
inputs         а
к "'в$
К
0         аГ
Ъ |
&__inference_dense_5_layer_call_fn_3376R4в1
*в'
%К"
inputs         а
к "К         аГд
A__inference_dense_6_layer_call_and_return_conditional_losses_3387_]^1в.
'в$
"К
inputs         аГ
к "&в#
К
0         Р
Ъ |
&__inference_dense_6_layer_call_fn_3396R]^1в.
'в$
"К
inputs         аГ
к "К         Рв
A__inference_dense_7_layer_call_and_return_conditional_losses_3407]de0в-
&в#
!К
inputs         Р
к "%в"
К
0         X
Ъ z
&__inference_dense_7_layer_call_fn_3416Pde0в-
&в#
!К
inputs         Р
к "К         Xв
?__inference_dense_layer_call_and_return_conditional_losses_3206_$%1в.
'в$
"К
inputs         аГ
к "&в#
К
0         Р
Ъ z
$__inference_dense_layer_call_fn_3215R$%1в.
'в$
"К
inputs         аГ
к "К         Рг
C__inference_dropout_2_layer_call_and_return_conditional_losses_3428\3в0
)в&
 К
inputs         X
p
к "%в"
К
0         X
Ъ г
C__inference_dropout_2_layer_call_and_return_conditional_losses_3433\3в0
)в&
 К
inputs         X
p 
к "%в"
К
0         X
Ъ {
(__inference_dropout_2_layer_call_fn_3438O3в0
)в&
 К
inputs         X
p
к "К         X{
(__inference_dropout_2_layer_call_fn_3443O3в0
)в&
 К
inputs         X
p 
к "К         Xг
A__inference_dropout_layer_call_and_return_conditional_losses_3307^4в1
*в'
!К
inputs         ╨A
p
к "&в#
К
0         ╨A
Ъ г
A__inference_dropout_layer_call_and_return_conditional_losses_3312^4в1
*в'
!К
inputs         ╨A
p 
к "&в#
К
0         ╨A
Ъ {
&__inference_dropout_layer_call_fn_3317Q4в1
*в'
!К
inputs         ╨A
p
к "К         ╨A{
&__inference_dropout_layer_call_fn_3322Q4в1
*в'
!К
inputs         ╨A
p 
к "К         ╨Aд
A__inference_flatten_layer_call_and_return_conditional_losses_3190_4в1
*в'
%К"
inputs         а
к "'в$
К
0         аГ
Ъ |
&__inference_flatten_layer_call_fn_3195R4в1
*в'
%К"
inputs         а
к "К         аГ╜
?__inference_model_layer_call_and_return_conditional_losses_2591z$%+,239:@AQR]^de=в:
3в0
&К#
input_1         ╨A
p

 
к "%в"
К
0         X
Ъ ╜
?__inference_model_layer_call_and_return_conditional_losses_2645z$%+,239:@AQR]^de=в:
3в0
&К#
input_1         ╨A
p 

 
к "%в"
К
0         X
Ъ ╝
?__inference_model_layer_call_and_return_conditional_losses_2985y$%+,239:@AQR]^de<в9
2в/
%К"
inputs         ╨A
p

 
к "%в"
К
0         X
Ъ ╝
?__inference_model_layer_call_and_return_conditional_losses_3077y$%+,239:@AQR]^de<в9
2в/
%К"
inputs         ╨A
p 

 
к "%в"
К
0         X
Ъ Х
$__inference_model_layer_call_fn_2741m$%+,239:@AQR]^de=в:
3в0
&К#
input_1         ╨A
p

 
к "К         XХ
$__inference_model_layer_call_fn_2836m$%+,239:@AQR]^de=в:
3в0
&К#
input_1         ╨A
p 

 
к "К         XФ
$__inference_model_layer_call_fn_3118l$%+,239:@AQR]^de<в9
2в/
%К"
inputs         ╨A
p

 
к "К         XФ
$__inference_model_layer_call_fn_3159l$%+,239:@AQR]^de<в9
2в/
%К"
inputs         ╨A
p 

 
к "К         Xг
A__inference_reshape_layer_call_and_return_conditional_losses_3335^0в-
&в#
!К
inputs         ╨A
к "*в'
 К
0         ╨A
Ъ {
&__inference_reshape_layer_call_fn_3340Q0в-
&в#
!К
inputs         ╨A
к "К         ╨A┤
"__inference_signature_wrapper_2879Н$%+,239:@AQR]^de@в=
в 
6к3
1
input_1&К#
input_1         ╨A"5к2
0
	dropout_2#К 
	dropout_2         X