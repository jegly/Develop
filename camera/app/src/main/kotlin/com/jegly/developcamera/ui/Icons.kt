package com.jegly.developcamera.ui

import androidx.compose.foundation.Canvas
import androidx.compose.foundation.layout.size
import androidx.compose.runtime.Composable
import androidx.compose.ui.Modifier
import androidx.compose.ui.geometry.CornerRadius
import androidx.compose.ui.geometry.Offset
import androidx.compose.ui.geometry.Size
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.Path
import androidx.compose.ui.graphics.drawscope.Stroke
import androidx.compose.ui.unit.Dp
import androidx.compose.ui.unit.dp
import kotlin.math.cos
import kotlin.math.sin

// Small hand-drawn icons, so the app needs no icon library. Line icons in the given colour.

@Composable
fun ApertureIcon(color: Color, size: Dp = 22.dp) = Canvas(Modifier.size(size)) {
    val r = this.size.minDimension / 2f
    val w = r * 0.14f
    drawCircle(color, r - w / 2, style = Stroke(w))
    // six blades: chords that leave a hexagonal opening in the middle
    for (i in 0 until 6) {
        val a = Math.toRadians(60.0 * i)
        val b = Math.toRadians(60.0 * i + 100.0)
        drawLine(
            color,
            Offset(center.x + (r * 0.92f * cos(a)).toFloat(), center.y + (r * 0.92f * sin(a)).toFloat()),
            Offset(center.x + (r * 0.32f * cos(b)).toFloat(), center.y + (r * 0.32f * sin(b)).toFloat()),
            strokeWidth = w,
        )
    }
}

@Composable
fun CanisterIcon(color: Color, size: Dp = 22.dp) = Canvas(Modifier.size(size)) {
    val s = this.size.minDimension
    val w = s * 0.08f
    drawRoundRect(color, Offset(s * 0.22f, s * 0.2f), Size(s * 0.46f, s * 0.68f), CornerRadius(s * 0.08f), style = Stroke(w))
    drawRect(color, Offset(s * 0.36f, s * 0.08f), Size(s * 0.18f, s * 0.12f))
    drawLine(color, Offset(s * 0.68f, s * 0.42f), Offset(s * 0.92f, s * 0.42f), strokeWidth = w)
    drawLine(color, Offset(s * 0.68f, s * 0.66f), Offset(s * 0.92f, s * 0.66f), strokeWidth = w)
}

@Composable
fun GearIcon(color: Color, size: Dp = 20.dp) = Canvas(Modifier.size(size)) {
    val r = this.size.minDimension / 2f
    val w = r * 0.18f
    drawCircle(color, r * 0.55f, style = Stroke(w))
    drawCircle(color, r * 0.18f)
    for (i in 0 until 8) {
        val a = Math.toRadians(45.0 * i)
        drawLine(
            color,
            Offset(center.x + (r * 0.62f * cos(a)).toFloat(), center.y + (r * 0.62f * sin(a)).toFloat()),
            Offset(center.x + (r * 0.95f * cos(a)).toFloat(), center.y + (r * 0.95f * sin(a)).toFloat()),
            strokeWidth = w * 1.3f,
        )
    }
}

/**
 * Front/back camera: a camera with two arrows going round it (jegly's pick "B", 2026-11-26). [hole] is the
 * colour behind the icon, for the lens.
 */
@Composable
fun FlipIcon(color: Color, hole: Color, size: Dp = 28.dp) = Canvas(Modifier.size(size)) {
    val u = this.size.minDimension * 0.86f   // the drawing's unit (the mockup's 24 dp icon)
    val c = center
    val w = u * 0.075f
    // The body, its top and the lens.
    val bw = u * 0.46f
    val bh = u * 0.32f
    val top = c.y - bh / 2 + u * 0.02f
    drawRoundRect(color, Offset(c.x - bw / 2, top), Size(bw, bh), CornerRadius(u * 0.05f))
    drawRect(color, Offset(c.x - bw * 0.18f, top - u * 0.05f), Size(bw * 0.36f, u * 0.07f))
    val lens = Offset(c.x, c.y + u * 0.02f)
    drawCircle(hole, u * 0.09f, lens)
    drawCircle(color, u * 0.025f, lens)
    // Two arrows round it, clockwise, each with its head at its end.
    val r = u * 0.47f
    val head = u * 0.11f
    for (start in floatArrayOf(205f, 25f)) {
        val sweep = 60f
        drawArc(color, start, sweep, false, Offset(c.x - r, c.y - r), Size(2 * r, 2 * r), style = Stroke(w))
        val t = Math.toRadians((start + sweep).toDouble())
        val ex = c.x + r * cos(t).toFloat()
        val ey = c.y + r * sin(t).toFloat()
        val tx = -sin(t).toFloat()   // the direction of travel at the end
        val ty = cos(t).toFloat()
        val nx = cos(t).toFloat()
        val ny = sin(t).toFloat()
        drawPath(
            Path().apply {
                moveTo(ex + tx * head, ey + ty * head)
                lineTo(ex - tx * head * 0.2f + nx * head * 0.8f, ey - ty * head * 0.2f + ny * head * 0.8f)
                lineTo(ex - tx * head * 0.2f - nx * head * 0.8f, ey - ty * head * 0.2f - ny * head * 0.8f)
                close()
            },
            color,
        )
    }
}
