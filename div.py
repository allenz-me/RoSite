
import re, pathlib

p = pathlib.Path('./content/posts')

md = list(p.rglob("*.md"))

# 将图片的相对路径修改为 Hugo 导出时的正确的相对路径。
changeFigure = lambda s: s.replace("../figures", "../../figures")


def changeFormula(s: str):
    """
    将公式里的 backslash 都换成 backslash，下划线里的 underscore 换成 backslash underscore
    代码块中的 underscore 保持不变！
    """
    f = lambda s: s.group().replace("\\", "\\\\").replace("_", "\\_").replace("\n\n", "\n").replace('*', "\\ast").replace('+', "\\+").replace('-', "\\-")
    s = re.sub(r"\$(.*?)\$", f, s)
    s = re.sub(r"\$\$((.|\n)+?)\$\$", f, s)
    return s
    
for i in md:
    text = i.open('r').read()
    text = changeFigure(text)
    text = changeFormula(text)
    with i.open('w') as f:
        f.write(text)
