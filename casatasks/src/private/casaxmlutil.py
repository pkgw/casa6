import os
import functools
import inspect
from xml.dom import minidom

import casatasks

# constants for generating converter method
__FUNCTION = 'override_args'
__ARGS = '_a'
__ARGS_DICT = '_d'
__ARGS_SUPPLIED = '_s'
__LOGLEVEL_IN_FUNCTION = 'INFO'

__DEBUG = False
if __DEBUG:
    from pprint import pprint


def xml_constraints_injector(func):
    """Decorator which loads constraints from a casatask XML file and apply them to the arguments of the decorated casatask.

    This method is designed as decorator for task methods. It executes as below:
    1. converts a constraints element of a CASA task XML to a Python code.
    2. evaluates the code, then a Python function is generated.
    3. executes the function and overrides task arguments to values defined by constraints tag.

    ex)
    a constraints tag of a CASA task XML:

    <constraints>
            <when param="timebin">
                <notequals type="string" value="">
                        <default param="timespan"><value type="string"/></default>
                </notequals>
            </when>
            <when param="fitmode">
                <equals value="list">
                        <default param="nfit"><value type="vector"><value>0</value></value></default>
                </equals>
                <equals value="auto">
                        <default param="thresh"><value>5.0</value></default>
                        <default param="avg_limit"><value>4</value></default>
                        <default param="minwidth"><value>4</value></default>
                        <default param="edge"><value type="vector"><value>0</value></value></default>
                </equals>
                <equals value="interact">
                        <default param="nfit"><value type="vector"><value>0</value></value></default>
                </equals>
            </when>
    </constraints>

    generated Python function code from the above XML:

    def override_args(_a, _d, _s):  # _a: position args based on *args
                                    # _d: dict[key: position name, val: corresponding position index of a key]
                                    #     to use to get a position index of args by position name
                                    # _s: boolean array, it is the same length as the position args,
                                    #     and the positions of user-supplied arguments are set to True
        if _d.get('timebin') is not None and _a[_d['timebin']] != '':
            if _d.get('timespan') is not None and _s[_d['timespan']] is False and _a[_d['timespan']] == "":
                _a[_d['timespan']] = ''
                casatasks.casalog.post("overrode argument: timespan -> ''", "INFO")
        if _d.get('fitmode') is not None and _a[_d['fitmode']] == 'list':
            if _d.get('nfit') is not None and _s[_d['nfit']] is False and _a[_d['nfit']] == "":
                _a[_d['nfit']] = [0]
                casatasks.casalog.post("overrode argument: nfit -> [0]", "INFO")
        if _d.get('fitmode') is not None and _a[_d['fitmode']] == 'auto':
            if _d.get('thresh') is not None and _s[_d['thresh']] is False and _a[_d['thresh']] == "":
                _a[_d['thresh']] = 5.0
                casatasks.casalog.post("overrode argument: thresh -> 5.0", "INFO")
            if _d.get('avg_limit') is not None and _s[_d['avg_limit']] is False and _a[_d['avg_limit']] == "":
                _a[_d['avg_limit']] = 4
                casatasks.casalog.post("overrode argument: avg_limit -> 4", "INFO")
            if _d.get('minwidth') is not None and _s[_d['minwidth']] is False and _a[_d['minwidth']] == "":
                _a[_d['minwidth']] = 4
                casatasks.casalog.post("overrode argument: minwidth -> 4", "INFO")
            if _d.get('edge') is not None and _s[_d['edge']] is False and _a[_d['edge']] == "":
                _a[_d['edge']] = [0]
                casatasks.casalog.post("overrode argument: edge -> [0]", "INFO")
        if _d.get('fitmode') is not None and _a[_d['fitmode']] == 'interact':
            if _d.get('nfit') is not None and _s[_d['nfit']] is False and _a[_d['nfit']] == "":
                _a[_d['nfit']] = [0]
                casatasks.casalog.post("overrode argument: nfit -> [0]", "INFO")

    Note: handling of <kwarg/> tag of task XML files
        Subparameters whose default value is the empty string '' - but where the empty string means in fact that 
        the real default value must be set to some non-empty string - require special care. One must be able to
        determine whether the empty string was user-supplied or not.
        To make this determination possible, the <kwarg/> tag must be set in the <param> tag definition of such
        parameters in the task XML file. This is currently the case only for parameter 'intent' of task sdcal,
        where intent='' means intent='all'. See sdcal.xml.


    Parameters
    ----------
    func : function
        The casatask function to be decorated

    Returns
    -------
    wrapper: function
        A decorated casatask satisfying the XML constraints
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        retval = None
        # Any errors are handled outside the task.
        # however, the implementation below is effectively
        # equivalent to handling it inside the task.
        funcname = func.__name__

        # load the function name and arguments which is wrapped the decorator
        # get an object reference to read informantion of argument
        func_ = func.__dict__.get('__wrapped__', func)

        is_recursive_load = False
        is_called_from_casatasks = False
        for frame_info in inspect.stack():
            if frame_info.function == func_.__name__:
                # when the task is called from the same task (ex: sdcal with two calmodes calls itself)
                is_recursive_load = True
            if frame_info.function == '__call__' and \
               frame_info.frame.f_locals['self'].__module__ == 'casatasks.' + func_.__name__:
                # check whether the function is called from casatasks or not.
                is_called_from_casatasks = True

        if is_recursive_load:
            casatasks.casalog.post('recursive task call', 'INFO')
            retval = func(*args, **kwargs)
        elif not is_called_from_casatasks:
            retval = func(*args, **kwargs)
        else:
            # generate the argument specification and the injector method from a task xml
            args_, args_position_dict, converter_function_string = __load_xml(funcname)

            # Note: the length of args is reduced by the length of kwargs.
            for i in range(len(args)):
                args_[i] = args[i]
            supplied_args_flags = [False] * len(args_position_dict)

            kwargs_ = dict()
            for k, v in kwargs.items():
                if args_position_dict.get(k) is not None:
                    args_[args_position_dict[k]] = v
                    supplied_args_flags[args_position_dict[k]] = True
                else:
                    kwargs_[k] = v

            if __DEBUG:
                print(converter_function_string)
                pprint(args_position_dict)
                pprint(args_)

            # override args by the converter generated from xml
            casatasks.casalog.post('loaded constraints from XML', 'DEBUG')
            _local = {}
            exec(converter_function_string, globals(), _local)
            _local[__FUNCTION](args_, args_position_dict, supplied_args_flags)

            # execute task
            retval = func(*args_, **kwargs_)

        return retval
    return wrapper
    

def __get_taskxmlfilepath(task):
    xmlpath = os.path.abspath(casatasks.__path__[0]) + '/__xml__'
    taskxmlfile = f'{xmlpath}/{task}.xml'
    if not os.path.isfile(taskxmlfile):
        raise ValueError
    if not os.access(taskxmlfile, os.R_OK):
        return PermissionError
    return taskxmlfile


def __load_xml(task):
    taskxml = __get_taskxmlfilepath(task)

    stmt = []
    dom = minidom.parse(taskxml)
    constraints = dom.getElementsByTagName('constraints')[0]
    for s in constraints.getElementsByTagName('when'):
        __handle_when(s, stmt)
    args = [__generate_default_value(param) for param in dom.getElementsByTagName('param')]
    args_position_dict = {param.getAttribute('name'): i for i, param in enumerate(dom.getElementsByTagName('param'))}
    return args, args_position_dict, __convert_stmt_to_pycode(stmt)


def __generate_default_value(param):
    type_ = param.getAttribute('type')
    value_, type_ = __handle_value(param.getElementsByTagName('value')[0], type_)
    if type_ == 'int':
        return int(value_)
    elif type_ == 'double':
        return float(value_)
    elif type_ == 'bool':
        return value_ == 'True'
    return value_


def __convert_stmt_to_pycode(stmt_list):
    ret = f'def {__FUNCTION}({__ARGS}, {__ARGS_DICT}, {__ARGS_SUPPLIED}):\n'
    if len(stmt_list) > 0:
        for [stmt, indent] in stmt_list:
            ret += __indent(indent) + stmt + '\n'
    else:
        ret += __indent(1) + 'pass\n'
    return ret


""" constants and methods for converting from XML tree to Python code """
__QUOTE = '\''
__OP_EQUALS = '=='
__OP_ASSIGN = '='
__OP_IS = 'is'
__OP_NOT_EQUAL = '!='
__OP_AND = 'and'
__OP_NOT = 'not'
__NONE = 'None'


def __handle_when(when, stmt):
    # <when>
    for node in when.childNodes:
        if node.nodeName == 'equals':
            __handle_equals_or_not_equals(when, node, stmt, __OP_EQUALS)
        elif node.nodeName == 'notequals':
            __handle_equals_or_not_equals(when, node, stmt, __OP_NOT_EQUAL)


def __handle_equals_or_not_equals(when, elem, stmt, operator):
    # <equals> or <notequals>
    indent_level = 1
    defaults = elem.getElementsByTagName('default')
    if len(defaults) > 0:
        stmt.append([__when(__get_param(when), operator, __get_value(elem)), indent_level])
        for default_ in defaults:
            left = __get_param(default_)
            right = default_.getElementsByTagName('value')[0]
            __handle_default(left, right, stmt, indent_level)


def __handle_default(left, right, stmt, indent_level):
    # <default>
    quote = ''
    right, type_ = __handle_value(right)
    if type_ == 'string' or type_ == 'record' or type_ == 'stringVec':
        quote = __QUOTE
    if type_.endswith('Vec') or type_ == 'vector':
        if isinstance(right, list):
            right = ','.join([f'{quote}{r}{quote}' for r in right])
        right = f'[{right}]'
    else:
        right = f'{quote}{right}{quote}'
    if_ = __if(
        __and(__can_get(__ARGS_DICT, left),
              __and(
                __is(__list(__ARGS_SUPPLIED, __dict(__ARGS_DICT, left)), False),
                __equals(__list(__ARGS, __dict(__ARGS_DICT, left)), '""'))
              )
        )
    stmt.append([if_, indent_level + 1])
    stmt.append([__assign(__list(__ARGS, __dict(__ARGS_DICT, left)), right), indent_level + 2])
    stmt.append([__casalog(left, right), indent_level + 2])


def __handle_value(_element, type_='int'):
    # <value>, it contains <value> tags or a value
    if _element.nodeName == 'value':
        if _element.hasAttribute('type'):
            type_ = _element.getAttribute('type')
        values = _element.getElementsByTagName('value')
        if len(values) > 0:
            return [__handle_value(v, type_)[0] for v in values], type_
        if _element.firstChild:
            return __handle_value(_element.firstChild, type_)
    elif hasattr(_element, 'data'):
        return _element.data, type_
    return '', type_


def __get_param(doc):
    return __get_attr(doc, 'param')


def __get_value(doc):
    if doc.hasAttribute('type') and doc.getAttribute('type') == 'vector':
        return __get_value(doc.firstChild)
    return __get_attr(doc, 'value')


def __get_attr(doc, param):
    s = doc.getAttribute(param)
    if s == '' or s:
        return s
    raise Exception('XML Parse Error')


def __when(left, operator, right):
    if ',' in right:
        right_ = ','.join(sorted([s.strip() for s in right.split(',')]))
        left_ = f"','.join(sorted([s.strip() for s in {__list(__ARGS, __dict(__ARGS_DICT, left))}.split(',')]))"
    else:
        right_ = right
        left_ = f'{__ARGS}[{__ARGS_DICT}[{__QUOTE}{left}{__QUOTE}]]'
    right_ = f'{__QUOTE}{right_}{__QUOTE}'

    return __if(__and(__can_get(__ARGS_DICT, left), __exp(left_, operator, right_)))


def __and(left, right):
    return __exp(left, __OP_AND, right)


def __assign(left, right):
    return __exp(left, __OP_ASSIGN, right)


def __is(left, right):
    return __exp(left, __OP_IS, right)


def __equals(left, right):
    return __exp(left, __OP_EQUALS, right)


def __exp(left, operator, right):
    return f'{left} {operator} {right}'


def __not(right):
    return f'{__OP_NOT} {right}'


def __if(exp):
    return f'if {exp}:'


def __get(val, operand, exp=None):
    if exp:
        return f'{val}.get({__QUOTE}{operand}{__QUOTE}, {exp})'
    return f'{val}.get({__QUOTE}{operand}{__QUOTE})'


def __can_get(val, operand):
    return __is(__get(val, operand), __not(__NONE))


def __casalog(left, right):
    return f'casatasks.casalog.post("overrode argument: {left} -> {right}", "{__LOGLEVEL_IN_FUNCTION}")'


def __dict(val, pos):
    return f'{val}[{__QUOTE}{pos}{__QUOTE}]'


def __list(val, pos):
    return f'{val}[{pos}]'


def __indent(level):
    return ' ' * 4 * level


if __name__ == '__main__':

    @xml_constraints_injector
    def sdcal(infile=None, calmode='tsys', fraction='10%', noff=-1,
              width=0.5, elongated=False, applytable='', interp='', spwmap={},
              outfile='', overwrite=False, field='', spw='', scan='', intent=''):
        print(calmode)
        print(fraction)
        print(intent)

    class _sdcal_py:

        def __call__(self, infile=None, calmode='tsys', fraction='10%', noff=-1,
                     width=0.5, elongated=False, applytable='', interp='', spwmap={},
                     outfile='', overwrite=False, field='', spw='', scan='', **kwargs):
            sdcal(infile, calmode, fraction, noff, width, elongated, applytable, interp, spwmap,
                  outfile, overwrite, field, spw, scan, **kwargs)

    # Note: on L124 "frame_info.frame.f_locals['self'].__module__" is '__main__' when below lines are executed,
    # so we cannot see the behavior of constraints in __main__ for now.
    # If you want to see it, replace the conditional expression to be True temporarily.
    x_sdcal = _sdcal_py()
    x_sdcal('test', calmode='otfraster,apply')
