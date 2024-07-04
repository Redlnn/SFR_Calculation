# -*- coding: utf-8 -*-  # noqa: UP009
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import cv2
import numpy as np

# from SFR import SFR

ROI_SIZE = 50
AOI_LENGTH = 30
AOI_WIDTH = 20


orig_img = cv2.imread((Path('edge_samples') / '0007.png').as_posix())

# is gray image?
try:
    gray_img = cv2.cvtColor(orig_img, cv2.COLOR_BGR2GRAY)
except Exception:
    gray_img = np.copy(orig_img)

# 图片尺寸
IMG_WIDTH = orig_img.shape[1]
IMG_HEIGHT = orig_img.shape[0]
diagonal_length = math.sqrt(IMG_WIDTH**2 + IMG_HEIGHT**2)

# 计算中心坐标
center_x = IMG_WIDTH // 2
center_y = IMG_HEIGHT // 2

cos_height = (IMG_HEIGHT**2 + diagonal_length**2 - IMG_WIDTH**2) / (2 * IMG_HEIGHT * diagonal_length)
cos_width = (IMG_WIDTH**2 + diagonal_length**2 - IMG_HEIGHT**2) / (2 * IMG_WIDTH * diagonal_length)


def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


@dataclass
class AOI:
    '''储存ROI中的T或S方向数据'''

    lt_pos: tuple[int, int]
    rt_pos: tuple[int, int]
    img: cv2.typing.MatLike


@dataclass
class ROI:
    '''储存各视场各方向数据'''

    img: cv2.typing.MatLike
    lt_pos: tuple[int, int]
    rt_pos: tuple[int, int]
    direction: Literal['C'] | Literal['UL'] | Literal['UR'] | Literal['BL'] | Literal['BR']
    field: float
    t_aoi: AOI | None = None
    s_aoi: AOI | None = None


# 计算中心视场需要计算MTF的区域
rois: list[ROI] = [
    ROI(
        img=gray_img[
            center_y - ROI_SIZE : center_y + ROI_SIZE,
            center_x - ROI_SIZE : center_x + ROI_SIZE,
        ],
        lt_pos=(center_x - ROI_SIZE, center_y - ROI_SIZE),
        rt_pos=(center_x + ROI_SIZE, center_y + ROI_SIZE),
        direction='C',
        field=0,
    )
]

# 获得每个视场每个方向需要计算MTF的区域
for field in (0.3, 0.5, 0.75, 0.9):
    # 计算偏移向量
    offset = field * (diagonal_length / 2)
    offset_x = int(offset * cos_width)
    offset_y = int(offset * cos_height)

    # 图像中间左上角坐标
    pos_ul = 'UL', (center_x - offset_x, center_y - offset_y)
    pos_ur = 'UR', (center_x + offset_x, center_y - offset_y)
    pos_bl = 'BL', (center_x - offset_x, center_y + offset_y)
    pos_br = 'BR', (center_x + offset_x, center_y + offset_y)

    for _direction, (cx, cy) in (pos_ul, pos_ur, pos_bl, pos_br):
        if (ltx := cx - ROI_SIZE) < 0:
            ltx = 0
        if (rtx := cx + ROI_SIZE) > IMG_WIDTH:
            rtx = IMG_WIDTH
        if (lty := cy - ROI_SIZE) < 0:
            lty = 0
        if (rty := cy + ROI_SIZE) > IMG_HEIGHT:
            rty = IMG_HEIGHT
        rois.append(
            ROI(
                img=gray_img[lty:rty, ltx:rtx],
                lt_pos=(ltx, lty),
                rt_pos=(rtx, rty),
                direction=_direction,
                field=field,
            )
        )

rect_size = 0
delta = 0

# 遍历ROI找角点，并以此找到AOI（单个TS方向）
for i in rois:

    '''
    shi-Tomas算法是对Harris角点检测算法的改进,一般会比Harris算法得到更好的角点.
    Harris算法的角点响应函数是将矩阵M的行列式值与M的迹相减,利用差值判断是否为角点,
    后来shi和Tomas提出改进的方法是,若矩阵M的两个特征值中较小的一个大于阈值,则认为
    他是角点,即:
        R=min(a1,a2)

    API:
        corners=cv.goodFeaturesToTrack(image, maxcorners, qualityLevel, minDistance)
        参数:
            image: 输入的灰度图像
            maxCorners: 获取角点数的数目
            qualityLevel: 该参数指出最低可接受的角点质量水平,在0~1之间
            minDistance: 角点之间的最小欧氏距离,避免得到相邻特征点
        返回:
            corners:搜索到的角点,在这里所有低于质量水平的角点被排除,然后把合格的角点按照质量排序,
            然后将质量较好的角点附近(小于最小欧氏距离)的角点删除,最后找到maxCorners个角点返回

    '''
    corners = cv2.goodFeaturesToTrack(
        i.img,
        4,
        0.3,
        5,
    )
    # TODO: 占用时间长
    # print(corners)
    # cv2.imshow('1', i.img)
    # cv2.waitKey()
    # exit()

    if i.direction == 'C':
        _xs = []
        _ys = []
        for _corner in corners:
            _x, _y = _corner.ravel()
            _xs.append(_x)
            _ys.append(_y)

        # test=========================
        for idx, _x in enumerate(_xs):
            cv2.circle(
                orig_img,
                (i.lt_pos[0] + int(_x), i.lt_pos[1] + int(_ys[idx])),
                2,
                (0, 0, 255) if idx == 0 else (255, 0, 0),
                2,
            )
        # test=========================

        distances: list[float] = []
        for idx in range(len(_xs)):
            for _ in range(idx + 1, len(_xs)):
                distances.append(distance(_xs[idx], _ys[idx], _xs[_], _ys[_]))

        # 棋盘大小
        rect_size = int(np.min(distances))

        AOI_LENGTH = int(rect_size / 3 * 2)
        AOI_WIDTH = int(rect_size / 5 * 2)

        _delta = []
        for idx in range(1, len(_xs)):
            _delta.append(_xs[idx] - _xs[0])
            _delta.append(_ys[idx] - _ys[0])

        delta = int(np.min(np.abs(_delta)))

        print(f'delta:{delta}')

    if corners is None:
        continue
    x, y = corners[0].ravel()
    x, y = int(x), int(y)

    # test=========================
    _xs = []
    _ys = []
    for _corner in corners:
        _x, _y = _corner.ravel()
        _xs.append(_x)
        _ys.append(_y)

    for idx, _x in enumerate(_xs):
        cv2.circle(
            orig_img,
            (i.lt_pos[0] + int(_x), i.lt_pos[1] + int(_ys[idx])),
            2,
            (0, 0, 255) if idx == 0 else (255, 0, 0),
            2,
        )
    # test=========================

    # 默认取所选点的上方和右方边缘
    # 防止ROI选出ROI识别区域
    # 竖直方向超出
    if (_lty := int(i.lt_pos[1] + y - (rect_size / 2) - AOI_LENGTH / 2)) > i.lt_pos[1]:
        if (_ltx := int(i.lt_pos[0] + x + (delta / 2) - AOI_WIDTH / 2)) < 0:
            _ltx = 0
        if (_rbx := int(i.lt_pos[0] + x + (delta / 2) + AOI_WIDTH / 2)) > IMG_WIDTH:
            _rbx = IMG_WIDTH
        if _lty < 0:
            _lty = 0
        if (_rby := int(i.lt_pos[1] + y - (rect_size / 2) + AOI_LENGTH / 2)) > IMG_HEIGHT:
            _rby = IMG_HEIGHT
        _lt_pos = (_ltx, _lty)
        _rb_pos = (_rbx, _rby)
    else:
        if (_ltx := int(i.lt_pos[0] + x - (delta / 2) - AOI_WIDTH / 2)) < 0:
            _ltx = 0
        if (_rbx := int(i.lt_pos[0] + x - (delta / 2) + AOI_WIDTH / 2)) > IMG_WIDTH:
            _rbx = IMG_WIDTH
        if (_lty := int(i.lt_pos[1] + y + (rect_size / 2) - AOI_LENGTH / 2)) < 0:
            _lty = 0
        if (_rby := int(i.lt_pos[1] + y + (rect_size / 2) + AOI_LENGTH / 2)) > IMG_HEIGHT:
            _rby = IMG_HEIGHT
        _lt_pos = (_ltx, _lty)
        _rb_pos = (_rbx, _rby)
    i.s_aoi = AOI(
        lt_pos=_lt_pos,
        rt_pos=_rb_pos,
        img=orig_img[_lt_pos[1] : _rb_pos[1], _lt_pos[0] : _rb_pos[0]],
    )

    # 水平方向超出
    if (_ltx := int(i.lt_pos[0] + x - (rect_size / 2) - AOI_LENGTH / 2)) > i.lt_pos[0]:
        if _ltx < 0:
            _ltx = 0
        if (_rbx := int(i.lt_pos[0] + x - (rect_size / 2) + AOI_LENGTH / 2)) > IMG_WIDTH:
            _rbx = IMG_WIDTH
        if (_lty := int(i.lt_pos[1] + y - (delta / 2) - AOI_WIDTH / 2)) < 0:
            _lty = 0
        if (_rby := int(i.lt_pos[1] + y - (delta / 2) + AOI_WIDTH / 2)) > IMG_HEIGHT:
            _rby = IMG_HEIGHT
        _lt_pos = (_ltx, _lty)
        _rb_pos = (_rbx, _rby)
    else:
        if (_ltx := int(i.lt_pos[0] + x + (rect_size / 2) - AOI_LENGTH / 2)) < 0:
            _ltx = 0
        if (_rbx := int(i.lt_pos[0] + x + (rect_size / 2) + AOI_LENGTH / 2)) > IMG_WIDTH:
            _rbx = IMG_WIDTH
        if (_lty := int(i.lt_pos[1] + y + (delta / 2) - AOI_WIDTH / 2)) < 0:
            _lty = 0
        if (_rby := int(i.lt_pos[1] + y + (delta / 2) + AOI_WIDTH / 2)) > IMG_HEIGHT:
            _rby = IMG_HEIGHT
        _lt_pos = (_ltx, _lty)
        _rb_pos = (_rbx, _rby)
    i.t_aoi = AOI(
        lt_pos=_lt_pos,
        rt_pos=_rb_pos,
        img=orig_img[_lt_pos[1] : _rb_pos[1], _lt_pos[0] : _rb_pos[0]],
    )

for i in rois:
    cv2.rectangle(orig_img, i.lt_pos, i.rt_pos, (0, 255, 255), 2)
    if i.t_aoi and i.t_aoi:
        cv2.rectangle(
            orig_img,
            i.t_aoi.lt_pos,
            i.t_aoi.rt_pos,
            (0, 255, 0),
            2,
        )
    if i.s_aoi and i.s_aoi:
        cv2.rectangle(
            orig_img,
            i.s_aoi.lt_pos,
            i.s_aoi.rt_pos,
            (0, 255, 0),
            2,
        )

cv2.namedWindow("result", 0)
cv2.resizeWindow("result", 1600, 1115)  # 设置窗口大小
cv2.imshow('result', orig_img)
# cv2.imwrite('test.png', rois[0].img)
cv2.waitKey()

# if TYPE_CHECKING:
#     assert rois[0].t_aoi and rois[0].s_aoi

# sfr = SFR(verbose=True)
# mtf, status = sfr.calc_sfr(cv2.cvtColor(rois[0].t_aoi.img, cv2.COLOR_BGR2GRAY))

# for i in mtf:
#     print(f'{i[0] / 0.0024}\t:{i[1]}')
# print(status)
